use std::borrow::Cow;
use std::cell::UnsafeCell;
use std::fs::File;
use std::io::{BufRead,BufReader,Write,BufWriter};
use std::mem::MaybeUninit;
use std::path::Path;

use ahash::AHashMap as HashMap;
use anyhow::{bail, Context, Result};
use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;

use crate::cli::Options;

pub fn check_dependencies(dependencies:&[&str]) -> Result<()> {
    for &executable in dependencies {
        if which::which(executable).is_err() {
            bail!("missing {executable} dependency, please check your system PATH");
        }
    }
    Ok(())
}

pub fn get_maxrss() -> f64 {
    
    let usage_self = unsafe {
        let mut usage = MaybeUninit::uninit();
        assert_eq!(libc::getrusage(libc::RUSAGE_SELF, usage.as_mut_ptr()), 0);
        usage.assume_init()
    };

    let usage_children = unsafe {
        let mut usage = MaybeUninit::uninit();
        assert_eq!(libc::getrusage(libc::RUSAGE_CHILDREN, usage.as_mut_ptr()), 0);
        usage.assume_init()
    };

    usage_self.ru_maxrss.max(usage_children.ru_maxrss) as f64 / (1024.0 * 1024.0)
}


pub fn get_file_reader(path: &Path) -> Box<dyn BufRead> {
    let file = match File::open(path) {
        Err(err) => panic!("Could not open file \"{}\": {}", path.display(), err),
        Ok(file) => file,
    };
    match path.extension() {
        Some(ext) if ext == "gz" => Box::new(BufReader::new(MultiGzDecoder::new(file))),
        _ => Box::new(BufReader::new(file)),
    }
}


pub fn get_file_writer(path: &Path) -> Box<dyn Write> {
    let file = match File::create(path) {
        Err(err) => panic!("Could not open/create file \"{}\": {}", path.display(), err),
        Ok(file) => file,
    };
    match path.extension() {
        Some(ext) if ext == "gz" => Box::new(
            BufWriter::new(GzEncoder::new(file, Compression::default()))
        ),
        _ => Box::new(BufWriter::new(file)),
    }
}


pub fn run_minimap2(reference: &Path, reads: &Path, bam: &Path, opts: &Options) -> Result<()> {

    use std::process::{Command, Stdio};

    let mm2_preset = match opts.mode {
        crate::cli::Mode::Hifi => "-xmap-hifi",
        crate::cli::Mode::Nano => "-xlr:hq",//"-xmap-ont",
    };
    
    let mm2_args = ["-t", &opts.nb_threads.to_string(), "-a", mm2_preset, "--no-long-join", "-r50", "-g2k", "-z200", reference.to_str().unwrap(), reads.to_str().unwrap()];
    let mut minimap2 = Command::new("minimap2")
        .args(mm2_args)
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn().context("cannot run minimap2")?;

    let samtools_args = ["sort", "--threads", &opts.nb_threads.to_string(), "-o", bam.to_str().unwrap()];
    let _samtools_sort = Command::new("samtools")
        .args(samtools_args)
        .stdin(Stdio::from(minimap2.stdout.take().unwrap()))
        .stderr(Stdio::null())
        .output()
        .expect("cannot run samtools sort");

    let status = minimap2.wait().context("minimap2 did not run successfully")?;
    if !status.success() {
        bail!("minimap2 exited with {status}");
    }

    // bam indexing

    let nb_threads = opts.nb_threads.max(4);
    let samtools_args = ["index", "-@", &nb_threads.to_string(), bam.to_str().unwrap()];
    let mut samtools_index = Command::new("samtools")
        .args(samtools_args)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .spawn().context("cannot run samtools index")?;
    
    let status = samtools_index.wait().context("samtools index did not run successfully")?;
    if !status.success() {
        bail!("samtools index exited with {status}");
    }

    Ok(())
}


pub fn insert_newlines(string:&str, every:usize) -> Cow<'_, str> {

    if string.len() <= every || every == 0 {
        return Cow::Borrowed(string)
    }

    let nb_newlines = (string.len()-1)/every;
    let mut res = String::with_capacity(string.len()+nb_newlines);
    for i in 0..nb_newlines {
        let beg = i * every;
        let end = (i+1) * every;
        res.push_str(&string[beg..end]);
        res.push('\n');
    }
    res.push_str(&string[nb_newlines*every..]);
    Cow::Owned(res)
}

// def insert_newlines(string:str, every:int=120) -> str:
//     return string if every <= 0 else '\n'.join(string[i:i+every] for i in range(0, len(string), every))


#[inline(always)]
pub fn counter_from_iter<T>(iter: impl Iterator<Item = T>) -> HashMap<T,usize> 
where
    T: std::cmp::Eq + std::hash::Hash
{    
    let mut counter = HashMap::new();
    for e in iter {
        counter.entry(e).and_modify(|e| *e += 1).or_insert(1);
    }
    counter
}


// from https://stackoverflow.com/a/65182786
#[derive(Copy, Clone)]
pub struct UnsafeSlice<'a, T> {
    slice: &'a [UnsafeCell<T>],
}
unsafe impl<T: Send + Sync> Send for UnsafeSlice<'_, T> {}
unsafe impl<T: Send + Sync> Sync for UnsafeSlice<'_, T> {}

impl<'a, T> UnsafeSlice<'a, T> {
    pub fn new(slice: &'a mut [T]) -> Self {
        let ptr = slice as *mut [T] as *const [UnsafeCell<T>];
        Self { slice: unsafe { &*ptr } }
    }
    
    /// # Safety
    /// It is UB if two threads write to the same index without synchronization.
    pub unsafe fn write(&self, i: usize, value: T) {
        let ptr = self.slice[i].get();
        unsafe { *ptr = value; }
    }
}