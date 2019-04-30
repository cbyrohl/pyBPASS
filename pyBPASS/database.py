import re as _re
import os as _os


class BPASSdatabase(object):

    def __init__(self, path, version):
        if not _os.path.isdir(path):
            raise ValueError(path + " is not a directory.")
        self.path = _os.path.abspath(path)
        self.version = version
        return

    def _zFromFilename(self, fname):
        zStr = _re.search("\\.z(.{3})\\.", fname).group(1)
        try:
            z = float("0."+zStr)
        except ValueError:
            z = float("1e-"+zStr[-1])
        return z
