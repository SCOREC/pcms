import pydiag.data.fieldlib
import pydiag.data.momlib
import pydiag.data.profile_data
import pydiag.data.nrgdata
import pydiag.data.srcmom_data


class RunDataFiles(object):
    """ Collection of the various output file objects of a GENE run

    The purpose of this class is to avoid rereading large binary files several times
    """

    def __init__(self, common):
        self.cm = common
        self.filedict = {"field": None, "nrg": None}
        for spec in common.specnames:
            self.filedict.update({"mom_{}".format(spec): None, "profile_{}".format(spec): None,
                                  "srcmom_{}".format(spec): None})

    def get_fileobject(self, name):
        if not self.filedict.get(name):
            if name == "field":
                self.filedict["field"] = pydiag.data.fieldlib.FieldFile(
                    "field{}".format(self.cm.fileextension), self.cm)
            elif name == "nrg":
                self.filedict["nrg"] = pydiag.data.nrgdata.NrgFile("nrg{}".format(self.cm.fileextension), self.cm)
                self.filedict["nrg"].generate_timeseries()
            elif name.startswith("mom_"):
                self.filedict[name] = pydiag.data.momlib.MomFile("{}{}".format(name, self.cm.fileextension),self.cm)
            elif name.startswith("profile_"):
                self.filedict[name] = pydiag.data.profile_data.ProfileFile("{}{}".format(name, self.cm.fileextension),self.cm)
            elif name.startswith("srcmom_"):
                self.filedict[name] = pydiag.data.srcmom_data.SrcmomFile("{}{}".format(name, self.cm.fileextension),self.cm)

            else:
                raise RuntimeError("Invalid file object identifier")
        return self.filedict[name]
