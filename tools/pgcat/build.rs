// With the `karan` feature, pgcat links karan's C library from the folder KARAN_LIB names (karan's `pixi run
// build-c` leaves libkaran.so and the runtime libraries it loads in build/c) and finds it there at run time.
fn main() {
    println!("cargo:rerun-if-env-changed=KARAN_LIB");
    if std::env::var_os("CARGO_FEATURE_KARAN").is_none() {
        return;
    }
    let dir = std::env::var("KARAN_LIB")
        .expect("the karan feature needs KARAN_LIB: the folder holding libkaran.so (karan's `pixi run build-c`, build/c)");
    println!("cargo:rustc-link-search=native={dir}");
    println!("cargo:rustc-link-lib=dylib=karan");
    println!("cargo:rustc-link-arg=-Wl,-rpath,{dir}");
}
