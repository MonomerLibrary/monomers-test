import gemmi
import sys

ccd_cif, list_cif = sys.argv[1:]

doc_ccd = gemmi.cif.read(ccd_cif)
doc_list = gemmi.cif.read(list_cif)

lib_ids = set(doc_list.find_block("comp_list").find_values("_chem_comp.id"))

ccd_ids = {}

for b in doc_ccd:
    mon = b.find_value("_chem_comp.id")
    st = b.find_value("_chem_comp.pdbx_release_status")
    ccd_ids.setdefault(st, []).append(mon)

ccd_ids = {x:set(ccd_ids[x]) for x in ccd_ids}

#for st in ccd_ids:
#    print(st, len(ccd_ids[st]), "in CCD")
#    print("lib-ccd", len(lib_ids - ccd_ids[st]))
#    print("ccd-lib", len(ccd_ids[st] - lib_ids))

print("REL in CCD:", len(ccd_ids["REL"]))
print(" included in monlib:", len(ccd_ids["REL"].intersection(lib_ids)))
print(" not included in monlib:", len(ccd_ids["REL"] - lib_ids))
print("OBS in CCD:", len(ccd_ids["OBS"]))
print(" included in monlib:", len(ccd_ids["OBS"].intersection(lib_ids)))
print("", " ".join(sorted(ccd_ids["OBS"].intersection(lib_ids))))

print("in monlib, but not in CCD:", len(lib_ids - (ccd_ids["REL"] | ccd_ids["OBS"])))
print(" ".join(sorted(lib_ids - (ccd_ids["REL"] | ccd_ids["OBS"]))))

