# hurp
scripts to create wrfchemi input files for WRF forecast
### Installation
WRFpy is installable via pip:
```
pip install git+https://github.com/ERA-URBAN/hurp
```

### Usage
hurp provides functionality depending on the used command-line switches:
```
usage: hurp [-h] [-c MY_CONFIG] --NEIpath NEIPATH --MACCpath MACCPATH
            --trafficpath TRAFFICPATH --wrfchemipath WRFCHEMIPATH
            [--wrfoutpath WRFOUTPATH] --outputdir OUTPUTDIR
            [--rundays RUNDAYS]
            datetime

HURP run script wrf chem Args that start with '--' (eg. --NEIpath) can also be
set in a config file ( or specified via -c). Config file syntax allows:
key=value, flag=true, stuff=[a,b,c] (for details, see syntax at
https://goo.gl/R74nmi). If an arg is specified in more than one place, then
commandline values override environment variables which override config file
values which override defaults.

positional arguments:
  datetime              Datetime string "%Y%m%d"

optional arguments:
  -h, --help            show this help message and exit
  -c MY_CONFIG, --my-config MY_CONFIG
                        config file path
  --NEIpath NEIPATH     Directory containing NEI files [env var: NEIpath]
  --MACCpath MACCPATH   Directory containing MACC files [env var: MACCpath]
  --trafficpath TRAFFICPATH
                        Directory containing traffic files [env var:
                        trafficpath]
  --wrfchemipath WRFCHEMIPATH
                        Directory containing wrfchemi files [env var:
                        wrfchemipath]
  --wrfoutpath WRFOUTPATH
                        Directory containing wrfout files with domain
                        information [env var: wrfoutpath]
  --outputdir OUTPUTDIR
                        Target output directory [env var: outputdir]
  --rundays RUNDAYS     number of WRF rundays

```
