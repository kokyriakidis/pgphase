#!/usr/bin/env python3
"""Build a diagnostic copy without changing production source or binary."""
from pathlib import Path
import hashlib, json, os, shlex, subprocess
root=Path('test_data/tmp_gap_fix41/trace-build')
root.mkdir(exist_ok=True)
(root/'compiler-tmp').mkdir(exist_ok=True)
original=Path('src/collect_phase.cpp').read_text()
s=original
start=s.index('std::optional<bool> complete_recovery_block_flip(')
end=s.index('\nstatic std::optional<SupportedAlleleEdge>',start)
helper=s[start:end]
helper=helper.replace('return std::nullopt;', '{ std::fprintf(stderr,"TRACE_FLIP_REJECT\\t%lld\\t%lld\\t%d\\n",(long long)upstream_phase_set,(long long)downstream_phase_set,__LINE__); return std::nullopt; }')
helper=helper.replace('    // Whole-block orientation needs evidence', '''    std::fprintf(stderr,"TRACE_VOTES\\t%lld\\t%lld\\t%d,%d,%d,%d\\t%d,%d\\t%d,%d,%d,%d\\t%d,%d;%d,%d\\t%d,%d\\n",(long long)upstream_phase_set,(long long)downstream_phase_set,counts[0][0],counts[0][1],counts[1][0],counts[1][1],snp_same,snp_cross,physical_counts[0][0],physical_counts[0][1],physical_counts[1][0],physical_counts[1][1],label_votes[0].first,label_votes[0].second,label_votes[1].first,label_votes[1].second,source_flip[0],source_flip[1]);
    // Whole-block orientation needs evidence''')
s=s[:start]+helper+s[end:]
s=s.replace('            consistent <= crossing[2]) return false;', '            consistent <= crossing[2]) { std::fprintf(stderr,"TRACE_PATH_CUT\\t%lld\\t%lld\\t%lld\\t%d,%d,%d\\n",(long long)chunk.candidates[static_cast<size_t>(first)].phase_set,(long long)chunk.candidates[static_cast<size_t>(sites[i])].key.sort_pos(),(long long)chunk.candidates[static_cast<size_t>(sites[i+1])].key.sort_pos(),crossing[0],crossing[1],crossing[2]); return false; }')

start=s.index('size_t stitch_complete_recovery_phase_blocks(')
end=s.index('\nbool msa_boundary_dropout_is_supported(',start)
f=s[start:end]
f=f.replace('    size_t joined = 0;', '''    for (size_t gi=0; gi<gauges.size(); ++gi) {
        const auto& g=gauges[gi];
        std::fprintf(stderr,"TRACE_MATRIX\\t%zu\\t%lld\\t%lld\\t%zu\\t%zu\\n",gi,(long long)g.beg,(long long)g.end,g.bam_sites.size(),g.bam_reads.size());
        for (size_t ci=0; ci<g.bam_sites.size(); ++ci) {
            const auto& v=g.bam_sites[ci];
            std::fprintf(stderr,"TRACE_SOURCE_SITE\\t%zu\\t%zu\\t%lld\\t%d\\t%d\\t%s\\t%lld\\t%d\\t%d\\t%d\\t%lld\\n",gi,ci,(long long)v.key.pos,(int)v.key.type,v.key.ref_len,v.key.alt.c_str(),(long long)v.phase_set,v.hap1_allele,v.hap2_allele,v.clean_snp,(long long)v.source_phase_set);
        }
        for (const auto& r:g.bam_reads) for (size_t oi=0; oi<r.observations.size(); ++oi)
            std::fprintf(stderr,"TRACE_SOURCE_OBS\\t%zu\\t%s\\t%d\\t%zu\\t%d\\t%d\\n",gi,r.qname.c_str(),r.mapq,r.observations[oi].first,r.observations[oi].second,oi<r.base_qualities.size()?r.base_qualities[oi]:0);
    }
    size_t joined = 0;''')
f=f.replace('        const auto imported =', '''        std::fprintf(stderr,"TRACE_SEAM\\t%lld\\t%lld\\t%lld\\t%lld\\t%lld\\t%lld\\n",(long long)window.beg,(long long)window.end,(long long)left_root,(long long)right_root,(long long)gauge.beg,(long long)gauge.end);
        const auto imported =''')
f=f.replace('        if (!valid || !has_bam_block || chain.size() < 2) continue;', '''        std::fprintf(stderr,"TRACE_CHAIN\\t%lld\\t%lld\\t%d\\t%d",(long long)window.beg,(long long)window.end,valid,has_bam_block);
        for (const auto ps:chain) std::fprintf(stderr,"\\t%lld",(long long)ps);
        std::fprintf(stderr,"\\n");
        if (!valid || !has_bam_block || chain.size() < 2) continue;''')
f=f.replace('        for (size_t i = 0; i < chain.size(); ++i)', '        const auto probe=complete_recovery_block_flip(chunk,gauge,window.left_phase_set,window.right_phase_set,initial_phase_set_candidates.at(window.left_phase_set),initial_phase_set_candidates.at(window.right_phase_set),std::min(opts.min_mapq,opts.recovery_min_mapq));\n        std::fprintf(stderr,"TRACE_PROBE_RESULT\\t%lld\\t%lld\\t%d\\n",(long long)window.beg,(long long)window.end,probe?(*probe?1:0):-1);\n        for (size_t i = 0; i < chain.size(); ++i)',1)
s=s[:start]+f+s[end:]
trace=root/'collect_phase_trace.cpp'; trace.write_text(s)
env=dict(os.environ,TMPDIR=str((root/'compiler-tmp').resolve()))
subprocess.run(['g++','-O3','-std=c++17','-Wall','-Wextra','-Isrc','-c',str(trace),'-o',str(root/'collect_phase_trace.o')],env=env,check=True)
line=next(l for l in Path('test_data/tmp_gap_next39/retained-build.log').read_text().splitlines() if l.startswith('g++') and '-o pgphase ' in l)
args=shlex.split(line)
args[args.index('pgphase')]=str(root/'pgphase-trace')
args[args.index('src/collect_phase.o')]=str(root/'collect_phase_trace.o')
subprocess.run(args,env=env,check=True)
Path('evaluations/2026-10-03-recovery-block-continuity/trace-build.json').write_text(json.dumps({'production_binary_sha256':hashlib.sha256(Path('pgphase').read_bytes()).hexdigest(),'production_phase_source_sha256':hashlib.sha256(original.encode()).hexdigest(),'trace_binary_sha256':hashlib.sha256((root/'pgphase-trace').read_bytes()).hexdigest(),'diagnostic_source':str(trace)},indent=2)+'\n')
