# GetOrganelle Database

### Installation from Github

By default, `get_organelle_config.py` will automatically access this repository to download and build the SeedDatabase and LabelDatabase of the latest version. e.g.

    get_organelle_config.py -a fungus_mt

Due to the unstable accessibility to Github in some regions, the `get_organelle_config.py` sometimes fails with connection error (e.g. timeout). In that case, try a couple of more times will simply work.

### Installation from local

If `get_organelle_config.py` still fails after many trials, supposing you want to install version `0.0.0` of `embplant_pt` and `embplant_mt`, download this repository and run `get_organelle_config.py` with the flag `--use-local`. Make your own database is feasible if you use the same directory structure, but not recommended. 

You can choose anyone of the following code blocks to install:
    
  1. Use `curl` to download the released database (ca. 80 MB):
    
    curl -L https://github.com/Kinggerm/GetOrganelleDep/releases/download/v0.0.0/v0.0.0.tar.gz | tar zx
    get_organelle_config.py -a embplant_pt,embplant_mt --use-local ./0.0.0
    
  2. Use svn to download part of the repository (ca. 80 MB):
  
    svn co https://github.com/Kinggerm/GetOrganelleDB/trunk/0.0.0
    get_organelle_config.py -a embplant_pt,embplant_mt --use-local ./0.0.0
    
  3. Use git clone to clone the entire repository (ca. 100 MB):
  
    git clone https://github.com/Kinggerm/GetOrganelleDB
    get_organelle_config.py -a embplant_pt,embplant_mt --use-local ./GetOrganelleDB/0.0.0
    
