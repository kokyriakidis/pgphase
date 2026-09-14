#!/usr/bin/env python3
"""Exercise genotype-preserving projection through both cache and pileup paths."""
import os,subprocess,tempfile,unittest
from pathlib import Path
import pysam
SCRIPT=Path(os.environ.get('PHASE_SCRIPT',Path(__file__).with_name('phase_vcf_from_hp.py')))


class Projection(unittest.TestCase):
    def test_multiallelic_genotype_and_exact_insertion(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp)
            header=pysam.VariantHeader();header.contigs.add('chr1',length=1000)
            header.formats.add('GT',1,'String','Genotype');header.formats.add('PS',1,'Integer','Phase set');header.add_sample('S')
            with pysam.VariantFile(str(d/'input.vcf'),'w',header=header) as v:
                records = [(100,('A','C','G'),(1,2)),(200,('A','ATT','AGG'),(1,2)),
                           (300,('A','C'),(0,1)),(400,('A','C'),(0,1)),
                           (500,('A','C'),(1,0)),(600,('A','C'),(1,0)),
                           (700,('A','C'),(1,0))]
                for pos,alleles,gt in records:
                    r=v.new_record(contig='chr1',start=pos,alleles=alleles);r.samples['S']['GT']=gt
                    if pos >= 500:
                        r.samples['S'].phased=True;r.samples['S']['PS']=501
                    v.write(r)
            pysam.tabix_compress(str(d/'input.vcf'),str(d/'input.vcf.gz'),force=True)
            pysam.tabix_index(str(d/'input.vcf.gz'),preset='vcf',force=True)
            with pysam.AlignmentFile(str(d/'reads.bam'),'wb',header={'HD':{'VN':'1.6'},'SQ':[{'SN':'chr1','LN':1000}]}) as b:
                for pos in (100,200,300):
                    for hp in (1,2):
                        for n in range(3):
                            r=pysam.AlignedSegment();r.query_name=f'{pos}_{hp}_{n}';r.reference_id=0;r.reference_start=pos;r.mapping_quality=60
                            if pos==200:
                                r.query_sequence='A'+('GG' if hp==1 else 'TT')+'ACGT';r.cigartuples=[(0,1),(1,2),(0,4)]
                            else:
                                r.query_sequence=({1:'G',2:'C'}[hp] if pos==100 else {1:'A',2:'C'}[hp])+'ACGT';r.cigartuples=[(0,5)]
                            r.set_tag('HP',hp);r.set_tag('PS',99);b.write(r)
            pysam.index(str(d/'reads.bam'))
            for cached in (False,True):
                output=d/f'out{cached}.vcf'
                cmd=['python3',str(SCRIPT),str(d/'reads.bam'),str(d/'input.vcf.gz'),str(output),'--region','chr1']
                if cached:cmd+=['--support-cache',str(d/'support.tsv'),'--rebuild-support-cache']
                subprocess.run(cmd,check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
                with pysam.VariantFile(str(output)) as v:
                    records=list(v)
                    self.assertEqual([r.samples['S']['GT'] for r in records],[(2,1),(2,1),(0,1),(0,1),(1,0),(1,0),(1,0)])
                    self.assertEqual([r.samples['S'].phased for r in records],[True,True,True,False,False,False,False])
                    self.assertEqual([r.samples['S']['PS'] for r in records],[99,99,99,None,None,None,None])

            native=d/'native.vcf'
            with pysam.VariantFile(str(native),'w',header=header) as v:
                r=v.new_record(contig='chr1',start=400,alleles=('A','C'))
                r.samples['S']['GT']=(1,0);r.samples['S'].phased=True;r.samples['S']['PS']=401;v.write(r)
            output=d/'fallback.vcf'
            cmd=['python3',str(SCRIPT),str(d/'reads.bam'),str(d/'input.vcf.gz'),str(output),
                 '--region','chr1','--native-vcf',str(native),'--retain-caller-phase-blocks',
                 '--fallback-site','401','--fallback-site','501']
            subprocess.run(cmd,check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
            with pysam.VariantFile(str(output)) as v:
                records=list(v)
                self.assertEqual(records[3].samples['S']['GT'],(1,0))
                self.assertTrue(records[3].samples['S'].phased)
                self.assertEqual(records[3].samples['S']['PS'],401)
                self.assertEqual(records[4].samples['S']['GT'],(1,0))
                self.assertTrue(records[4].samples['S'].phased)
                self.assertEqual(records[4].samples['S']['PS'],1000000501)
                self.assertFalse(records[5].samples['S'].phased)


if __name__=='__main__':unittest.main()
