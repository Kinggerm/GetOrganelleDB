# GetOrganelle Database

All versions of default databases of [GetOrganelle](https://github.com/Kinggerm/GetOrganelle).

### Default directory

The initialized database will be located at `~/.GetOrganelle` by default. It can be changed by using the command line parameter `--config-dir` for a single run, or by using the shell environment variable `GETORG_PATH` for the entire running environment.

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

Due to the unstable accessibility to Github in some regions, the `get_organelle_config.py` sometimes fails with connection error (e.g. timeout). In most cases, try the above command for more times will simply work.

### Option 2 Initialization from local files

If `Initialization from Github` still fails after many trials, download this repository and run `get_organelle_config.py` with the flag `--use-local`. Making your own database is feasible if you use the same directory structure, but not recommended. 

Supposing you want to install version `0.0.0` of `embplant_pt` and `embplant_mt`, you can choose any one of the following code blocks to install:
    
  1. Use `curl` to download the released compressed file (ca. 20 MB -> 80 MB):
    
    curl -L https://github.com/Kinggerm/GetOrganelleDep/releases/download/v0.0.1/v0.0.1.tar.gz | tar zx
    get_organelle_config.py -a embplant_pt,embplant_mt --use-local ./0.0.1
    
  2. Use svn to download part of the repository (ca. 80 MB):
  
    svn co https://github.com/Kinggerm/GetOrganelleDB/trunk/0.0.1
    get_organelle_config.py -a embplant_pt,embplant_mt --use-local ./0.0.1
    
  3. Use git clone to clone the entire repository (ca. 200 MB):
  
    git clone https://github.com/Kinggerm/GetOrganelleDB
    get_organelle_config.py -a embplant_pt,embplant_mt --use-local ./GetOrganelleDB/0.0.1
    

### Updates

* **0.0.1** fungus_nr added

* **0.0.0** Initial version
