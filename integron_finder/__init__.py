import sys

__version__ = '$VERSION'

class IntegronError(Exception):
    pass


def get_version_message():
    # if week keep '$VERSION' as is
    # the setup.py will replace it by the value set in setup
    # so the test become True even if integron_finder is installed using setup.py
    if __version__ == '$' + 'VERSION':
        version = "NOT packaged, it should be development version"
    else:
        version = __version__
    version_text = """integron_finder version {0}
Python {1}

 - open-source GPLv3,
 - Jean Cury, Bertrand Neron, Eduardo Rocha,
 - citation:

 Identification and analysis of integrons and cassette arrays in bacterial genomes
 Jean Cury; Thomas Jove; Marie Touchon; Bertrand Neron; Eduardo PC Rocha
 Nucleic Acids Research 2016; doi: 10.1093/nar/gkw319
 """.format(version, sys.version)
    return version_text