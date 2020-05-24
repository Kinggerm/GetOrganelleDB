# GetOrganelle Database

By default, `get_organelle_config.py` will automatically access this repository to download and build the SeedDatabase and LabelDatabase of the latest version.

    get_organelle_config.py -a fungus_mt

Due to the unstable accessibility to Github in some regions, the `get_organelle_config.py` sometimes fails with connection error (e.g. timeout). In that case, try a couple of more times will simply work.

If `get_organelle_config.py` still fails after many trials, supposing you want to install version `0.0.0` of `embplant_pt` and `embplant_mt`, download this repository and run `get_organelle_config.py` with the flag `--use-local`:
    
    version=0.0.0
    git clone https://github.com/Kinggerm/GetOrganelleDB
    get_organelle_config.py -a embplant_pt,embplant_mt --use-local ./GetOrganelleDB/$version
    
    # or download part of the repository by
    # svn co https://github.com/Kinggerm/GetOrganelleDB/trunk/$version
    # get_organelle_config.py -a embplant_pt,embplant_mt --use-local ./$version
    
