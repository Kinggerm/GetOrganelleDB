# GetOrganelle Database

All versions of default databases of [GetOrganelle](https://github.com/Kinggerm/GetOrganelle)

### Default directory

By default, the initialized database will be located at `~/.GetOrganelle`. It can be changed, by using the command line parameter `--config-dir` for a single run, or by using the shell environment variable `GETORG_PATH` for the entire running environment.

For example, one may change the default database directory into `/home/shared/.GetOrganelle` by adding 

> GETORG_PATH=/home/shared/.GetOrganelle<br>
> export GETORG_PATH

to `/etc/profile` for **system-wide** usage, 
or to `~/.bashrc` for **user-wide** usage in **Ubuntu Desktop**, 
or to `~/.bash_profile` for **user-wide** usage in **bash**, 
or to `~/.zshrc` for **user-wide** usage in **Zsh**

and restarting the shell before initialization.

### Option 1 Initialization from Github

By default, `get_organelle_config.py` will automatically access this repository to download and build the SeedDatabase and LabelDatabase of the latest version. e.g.

    get_organelle_config.py -a fungus_mt

Due to the unstable accessibility to Github in some regions, the `get_organelle_config.py` sometimes fails with connection error (e.g., timeout, sha256_unmatch). However, trying the above command more times will simply work in most cases.

### Option 2 Initialization from local files

If `Initialization from Github` still fails after many trials, download this repository and run `get_organelle_config.py` with the flag `--use-local`. Making your own database is feasible if you use the same directory structure but not recommended. 

Supposing you want to install version `0.0.1` of `embplant_pt` and `embplant_mt`, you can choose any one of the following code blocks to install:
    
  1. Use `curl` to download the released compressed file (ca. 20 MB -> 80 MB):
    
    curl -L https://github.com/Kinggerm/GetOrganelleDep/releases/download/v0.0.1/v0.0.1.tar.gz | tar zx
    get_organelle_config.py -a embplant_pt,embplant_mt --use-local ./0.0.1
    
  2. Use svn to download part of the repository (ca. 80 MB):
  
    svn co https://github.com/Kinggerm/GetOrganelleDB/trunk/0.0.1
    get_organelle_config.py -a embplant_pt,embplant_mt --use-local ./0.0.1
    
  3. Use git clone to clone the entire repository (ca. 200 MB):
  
    git clone https://github.com/Kinggerm/GetOrganelleDB
    get_organelle_config.py -a embplant_pt,embplant_mt --use-local ./GetOrganelleDB/0.0.1
    

## Updates

* **0.0.1.minima** A minimal subset of 0.0.1 for test only (GetOrganelle 1.7.6+ required).

* **0.0.1** fungus_nr added

* **0.0.0** Initial version



## How to contribute

### Welcome to the community!

For generating previous databases, I downloaded all available **complete** plastomes/mitogenomes or **complete** nr region from the Genbank, semi-manually cleaned them, masked the simple repeats. However, using all sequences as the default database takes too much space for most users. So I made a customized script, `scripts/generate_bowtie2_seed.py`, to balance the taxa coverage and the size. The basic idea is randomly picking a set of sequences that can represent/recruit all the candidate sequences given a certain gapping threshold (default: 2000 bp). Then I used the sequence as the seed database and extracted the non-tRNA genes as the label database. Of course, each step should be accompanied by error-proof checking.

### Guidelines

1. The sequence header of the seed database should include the Genbank accession numbers.
2. The sequence header of the label database should be in the form of `>name type - sequence_info_without_space`, e.g. `>rbcL gene - Amborella_AJ506156_2`. The sequence info should at least include the Genbank accession number.
3. For the seed database, single- and oligo-nucleotide repeats should be masked using `N`s to alleviate the computational burden from calling irrelevant reads.
4. After using `scripts/generate_bowtie2_seed.py`, manually check and replace the representative sequences.
5. For the label database, remember to unify the gene names (e.g., upper/lower cases, abbreviations) and exclude the short genes such as tRNAs.
6. Test the database with a range of WGS from different taxa.

### Uploading
Follow the following steps to add the compiled database to the community:
1. fork the [GetOrganelleDB repo](https://github.com/Kinggerm/GetOrganelleDB).
2. duplicate the subdirectory of the latest version and rename it to a newer one.
   > e.g. if the latest is 0.0.1, rename the duplicate as 0.0.2
3. add the compiled seed/label databased inside the new subdirectory. 
4. send a pull request.

