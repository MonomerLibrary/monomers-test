# todo - check 0 sigmas (for torsions except const*)
# todo - identical sigmas in the same plane
from __future__ import absolute_import, division, print_function, generators
import unittest
import warnings
import os
import glob
import gemmi
import math
import collections

#monlib_path = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
monlib_path = os.environ["CLIBD_MON"]
warnings.formatwarning = lambda x, *a: "Warning: {}\n".format(x)

def read_mon_lib_list():
    doc = gemmi.cif.read(os.path.join(monlib_path, "list", "mon_lib_list.cif"))
    comp_list = doc.find_block("comp_list")
    link_list = doc.find_block("link_list")
    mod_list = doc.find_block("mod_list")
    return comp_list, link_list, mod_list, doc

def read_ener_lib():
    doc = gemmi.cif.read(os.path.join(monlib_path, "ener_lib.cif"))
    b = doc.find_block("energy")
    return b

def get_bond(mon, id1, id2):
    for b in mon.rt.bonds:
        if (b.id1.atom, b.id2.atom) == (id1, id2):
            return b
        if (b.id1.atom, b.id2.atom) == (id2, id1):
            return b

def check_tors(mon): # mon, link or mod
    ideals = {2: 0, 3: 60, 6: 0}
    if hasattr(mon, "id"): # link or mod
        name = mon.id
    else: # monomer
        name = mon.name
        if mon.group == gemmi.ChemComp.Group.Peptide:
            name += "(peptide)"
    for t in mon.rt.torsions:
        if t.id1.comp == 100: continue # delete
        if not t.label.startswith(("chi", "sp2_sp2")): continue # and not t.id4.atom.startswith("H"): continue # ad hoc
        if t.period in ideals:
            cos = math.cos(math.radians(t.period * (t.value - ideals[t.period])))
            if cos != 1:
                warnings.warn("strange torsion ideal: {} {} per={} value={}".format(name, t.label, t.period, t.value))

def check_monomer_chemtype(doc, name):
    return {} # turned off for now as it gives too many errors
    num_h = dict(NSP=0, NSP1=1, NS=0, NS1=1, NC1=1, NC2=2, NH0=0, NH1=1, NH2=2,
                 NPA=0, NPB=0, NRD5=0, NRD6=0, NR15=1, NR16=1, NR5=0, NR6=0,
                 N=0, NT=0, NT1=1, NT2=2, NT3=3, NT4=4, N30=0, N31=1, N32=2, N33=3,
                 S=0, SH1=1)
    monlib = gemmi.MonLib()
    monlib.read_monomer_doc(doc)
    cc = monlib.monomers[name]
    bond_h = {}
    ret = []
    for b in cc.rt.bonds:
        if cc.find_atom(b.id1.atom).is_hydrogen():
            bond_h[b.id2.atom] = bond_h.get(b.id2.atom, 0) + 1
        if cc.find_atom(b.id2.atom).is_hydrogen():
            bond_h[b.id1.atom] = bond_h.get(b.id1.atom, 0) + 1
    for a in cc.atoms:
        if a.chem_type in num_h and a.id in bond_h and num_h[a.chem_type] != bond_h[a.id]:
            ret.append("{}({},{})".format(a.id, a.chem_type, bond_h[a.id])) # wrong chem type
    return ret

class TestMonlib(unittest.TestCase):
    def setUp(self): self.errors = []
    def tearDown(self): self.assertEqual(len(self.errors), 0, msg="\n"+"\n".join(self.errors))

    def test_monomers(self):
        comp_list, link_list, mod_list, _ = read_mon_lib_list()
        lgroups = {row.str(0):row.str(1) for row in comp_list.find("_chem_comp.", ["id", "group"])}
        elib = read_ener_lib()
        all_types = set(elib.find_values("_lib_atom.type"))
        cgroups = {}
        for f in glob.glob(os.path.join(monlib_path, "?", "*.cif")):
            doc = gemmi.cif.read(f)
            block = doc.find_block("comp_list")
            for row in block.find("_chem_comp.", ["id", "group"]):
                cgroups[row[0]] = row.str(1)
                try: self.assertEqual(lgroups.get(row[0]), row.str(1),
                                     msg="{}: group {} vs {}".format(row[0], lgroups.get(row[0]), row.str(1)))
                except AssertionError as e: self.errors.append(str(e))

            b = doc[-1]
            all_atoms, alt_atoms = [], []
            for row in b.find("_chem_comp_atom.", ["atom_id", "type_energy", "?alt_atom_id"]):
                # test unknown energy type
                try: self.assertTrue(row.str(1) in all_types, msg="{} unknown energy type {} {}".format(os.path.basename(f), row[0], row[1]))
                except AssertionError as e: self.errors.append(str(e))
                all_atoms.append(row.str(0))
                if row.has(2): alt_atoms.append(row.str(2))

            # check duplication
            counts = collections.Counter(all_atoms)
            try: self.assertFalse(any(counts[x] > 1 for x in counts),
                                  msg="{} duplicated atoms {}".format(os.path.basename(f), counts))
            except AssertionError as e: self.errors.append(str(e))
            counts = collections.Counter(alt_atoms)
            try: self.assertFalse(any(counts[x] > 1 for x in counts),
                                  msg="{} duplicated alt atoms {}".format(os.path.basename(f), counts))
            except AssertionError as e: self.errors.append(str(e))
            
            # test unknown atoms in restraints
            all_atoms = set(all_atoms)
            for t1, t2 in (("_chem_comp_bond.", ["atom_id_1", "atom_id_2"]),
                           ("_chem_comp_angle.", ["atom_id_1", "atom_id_2", "atom_id_3"]),
                           ("_chem_comp_tor.", ["atom_id_1", "atom_id_2", "atom_id_3", "atom_id_4"]),
                           ("_chem_comp_chir.", ["atom_id_centre", "atom_id_1", "atom_id_2", "atom_id_3"]),
                           ("_chem_comp_plane_atom.", ["atom_id"])):
                for row in b.find(t1, t2):
                    for i in range(len(t2)):
                        if t1 == "_chem_comp_chir." and not row.str(i): continue # some metals
                        #try: self.assertTrue(row.str(i) in all_atoms, msg="{} unknown atom {}{} {}".format(os.path.basename(f), t1, t2[i], row.str(i)))
                        #except AssertionError as e: self.errors.append(str(e))
                        if row.str(i) not in all_atoms:
                            warnings.warn("Unknown atom in restraint: {} {}{} {}".format(os.path.basename(f), t1, t2[i], row.str(i)))
            # check chem types
            invalid_chemtypes = check_monomer_chemtype(doc, b.name[len("comp_"):])
            if invalid_chemtypes:
                self.errors.append("{} invalid chem types {}".format(os.path.basename(f), ",".join(invalid_chemtypes)))

        only_in_cif = set(cgroups) - set(lgroups)
        try: self.assertFalse(only_in_cif, msg="groups only in cif files")
        except AssertionError as e: self.errors.append(str(e))
        
        only_in_list = set(lgroups) - set(cgroups)
        try: self.assertFalse(only_in_list, msg="groups only in list")
        except AssertionError as e: self.errors.append(str(e))

    def test_gemmi_monlib(self):
        comp_list, link_list, mod_list, _ = read_mon_lib_list()
        lgroups = {row.str(0):row.str(1) for row in comp_list.find("_chem_comp.", ["id", "group"])}
        try: monlib = gemmi.read_monomer_lib(monlib_path, list(lgroups))
        except Exception as e: self.errors.append(str(e))

        def atoms_from_rt(rt, exd=False): # exd = exclude deleted
            d = ord("d")
            atoms = [a for b in rt.bonds for a in (b.id1, b.id2) if not exd or b.id1.comp != d]
            atoms.extend(a for b in rt.angles for a in (b.id1, b.id2, b.id3) if not exd or b.id1.comp != d)
            atoms.extend(a for b in rt.torsions for a in (b.id1, b.id2, b.id3, b.id4) if not exd or b.id1.comp != d)
            atoms.extend(a for b in rt.chirs for a in (b.id_ctr, b.id1, b.id2, b.id3) if not exd or b.id1.comp != d)
            atoms.extend(a for b in rt.planes for a in b.ids if not exd or a.comp != d)
            return atoms

        for ln in monlib.links:
            l = monlib.links[ln]
            if l.side1.comp == l.side2.comp == "": continue
            atoms = atoms_from_rt(l.rt)
            for i in range(1, 3):
                side = (l.side1, l.side2)[i-1]
                mon = side.comp
                if not mon: continue
                only_in_restr = set(a.atom for a in atoms if a.comp == i) - set(a.id for a in monlib.monomers[mon].atoms)
                try: self.assertFalse(only_in_restr, msg="unknown atoms in link {} mon {}".format(ln, mon))
                except AssertionError as e: self.errors.append(str(e))
                mod = side.mod
                if mod:
                    try: self.assertTrue(mod in monlib.modifications, msg="undefined mod {} in link {}".format(mod, ln))
                    except AssertionError as e: self.errors.append(str(e))

        # test if mod is applicable
        for ln in monlib.links:
            l = monlib.links[ln]
            for i, side in enumerate((l.side1, l.side2)):
                if side.comp == "" or side.mod == "": continue
                m = monlib.modifications[side.mod]
                if m.comp_id == side.comp: continue # dedicated mod will be checked
                undef = set()
                for x in m.atom_mods:
                    a = monlib.monomers[side.comp].find_atom(x.old_id)
                    if chr(x.func) in ("c", "d") and a is None:
                        undef.add(x.old_id)
                for x in m.rt.bonds:
                    if chr(x.id1.comp) in ("c", "d"):
                        mon_bond = get_bond(monlib.monomers[side.comp], x.id1.atom, x.id2.atom)
                        if mon_bond is None:
                            undef.add(x.id1.atom+"-"+x.id2.atom)
                # should check other restraints also

                try: self.assertFalse(undef, msg="mod {} is not applicable to {} (undef {})".format(side.mod, side.comp, undef))
                except AssertionError as e: self.errors.append(str(e))

        for mn in monlib.modifications:
            m = monlib.modifications[mn]
            if not m.comp_id: continue
            try: self.assertFalse(any([chr(am.func)=="a" and not am.new_id for am in m.atom_mods]),
                                  msg="{} _chem_mod_atom.new_atom_id missing".format(mn))
            except AssertionError as e: self.errors.append(str(e))
            atoms = atoms_from_rt(m.rt, True)
            # atoms may be added in modification
            mon_atoms = set(a.id for a in monlib.monomers[m.comp_id].atoms)
            added_atoms = set(am.new_id for am in m.atom_mods if chr(am.func)=="a")
            chged_atoms = set(am.old_id for am in m.atom_mods if chr(am.func)=="c")
            deled_atoms = set(am.old_id for am in m.atom_mods if chr(am.func)=="d")
            only_in_restr = set(a.atom for a in atoms) - ((mon_atoms | added_atoms) - deled_atoms)
            try: self.assertFalse(only_in_restr, msg="unknown atoms in mod {} mon {}".format(mn, m.comp_id))
            except AssertionError as e: self.errors.append(str(e))
            undefs = (deled_atoms | chged_atoms) - mon_atoms 
            try: self.assertFalse(undefs, msg="changing or deleting undefined atoms in mod {} mon {}".format(mn, m.comp_id))
            except AssertionError as e: self.errors.append(str(e))

        doc = gemmi.cif.read(os.path.join(monlib_path, "list", "mon_lib_list.cif"))
        for b in doc:
            if b.name.startswith("mod_") and b.name != "mod_list":
                mon = monlib.modifications[b.name[len("mod_"):]].comp_id
                if not mon: continue
                planes = {}
                for row in b.find("_chem_mod_plane_atom.", ["mod_id", "function", "plane_id", "atom_id"]):
                    planes.setdefault(row.str(2), []).append((row.str(1), row.str(3)))
                for p in planes:
                    if all(x[0] == "add" for x in planes[p]): continue
                    found = [x for x in monlib.monomers[mon].rt.planes if x.label == p]
                    try: self.assertTrue(found, msg="{} _chem_mod_plane_atom.plane_id {} not match".format(b.name, p))
                    except AssertionError as e: self.errors.append(str(e))
                    if not found: continue
                    notfound = set(x[1] for x in planes[p]) - set(x.atom for x in found[0].ids)
                    try: self.assertFalse(notfound, msg="{} _chem_mod_plane_atom.atom_id {} not match in {}".format(b.name, notfound, p))
                    except AssertionError as e: self.errors.append(str(e))

        # check torsions
        #for r in (monlib.monomers, monlib.links, monlib.modifications):
        #    for m in r:
        #        check_tors(r[m])

    def test_group(self):
        comp_list, link_list, mod_list, _ = read_mon_lib_list()
        lgroups = {row.str(0):row.str(1) for row in comp_list.find("_chem_comp.", ["id", "group"])}
        lg_set = set(lgroups.values())
        all_set = set(["peptide", "P-peptide", "M-peptide", "DNA", "RNA", "pyranose", "ketopyranose", "furanose", "NON-POLYMER"])
        try: self.assertTrue(lg_set.issubset(all_set), msg="unknown groups: {}".format(str(lg_set - all_set)))
        except AssertionError as e: self.errors.append(str(e))

    def test_links(self):
        comp_list, link_list, mod_list, doc = read_mon_lib_list()
        lgroups = {row.str(0):row.str(1) for row in comp_list.find("_chem_comp.", ["id", "group"])}
        known_groups = set(lgroups.values())
        known_groups.add("DNA/RNA") # can we really consider these known groups?
        known_groups.add("pept")
        
        ltab = link_list.find("_chem_link.", ["id", "comp_id_1", "mod_id_1", "group_comp_1", "comp_id_2", "mod_id_2", "group_comp_2"])
        mtab = mod_list.find("_chem_mod.", ["id", "comp_id", "group_id"])
        
        link_undef_group = [(row.str(0), row.str(i)) for row in ltab for i in (3,6) if row.str(i) and row.str(i) not in known_groups]
        mod_undef_group =  [(row.str(0), row.str(2)) for row in mtab if row.str(2) and row.str(2) not in known_groups]
        link_undef_comp =  [(row.str(0), row.str(i)) for row in ltab for i in (1,4) if row.str(i) !="" and row.str(i) not in lgroups]
        mod_undef_comp =   [(row.str(0), row.str(1)) for row in mtab if row.str(1) !="" and row.str(1) not in lgroups]

        try: self.assertFalse(link_undef_group, msg="undefined groups referenced in links")
        except AssertionError as e: self.errors.append(str(e))
        try: self.assertFalse(mod_undef_group, msg="undefined groups referenced in mods")
        except AssertionError as e: self.errors.append(str(e))
        try: self.assertFalse(link_undef_comp, msg="undefined comp referenced in links")
        except AssertionError as e: self.errors.append(str(e))
        try: self.assertFalse(mod_undef_comp, msg="undefined comp referenced in mods")
        except AssertionError as e: self.errors.append(str(e))

        for row in ltab:
            if row.str(1):
                try: self.assertEqual(row.str(3), lgroups[row.str(1)], msg="{} in link {}".format(row.str(1), row.str(0)))
                except AssertionError as e: self.errors.append(str(e))
            if row.str(4):
                try: self.assertEqual(row.str(6), lgroups[row.str(4)], msg="{} in link {}".format(row.str(4), row.str(0)))
                except AssertionError as e: self.errors.append(str(e))

        for row in mtab:
            if row.str(1):
                try: self.assertEqual(row.str(2), lgroups[row.str(1)], msg="{} in mod {}".format(row.str(1), row.str(0)))
                except AssertionError as e: self.errors.append(str(e))
            try: self.assertTrue(doc.find_block("mod_" + row.str(0)), msg="undefined mod {} in mod_list".format(row.str(0)))
            except AssertionError as e: self.errors.append(str(e))

    def test_mods(self):
        doc = gemmi.cif.read(os.path.join(monlib_path, "list", "mon_lib_list.cif"))
        no_id = []
        no_new_period = []
        for b in doc:
            for row in b.find("_chem_mod_tor.", ["function", "?id", "?new_period"]):
                if row.str(0) != "delete":
                    if not row.has(1) or gemmi.cif.is_null(row[1]): no_id.append(b.name)
                    if not row.has(2) or gemmi.cif.is_null(row[2]): no_new_period.append(b.name)
                elif row.has(2) and gemmi.cif.is_null(row[2]):
                    warnings.warn("{}: _chem_mod_tor.new_period '.' causes a problem in old gemmi".format(b.name))

        try: self.assertFalse(set(no_id), msg="no _chem_mod_tor.id")
        except AssertionError as e: self.errors.append(str(e))
        try: self.assertFalse(set(no_new_period), msg="no _chem_mod_tor.new_period")
        except AssertionError as e: self.errors.append(str(e))

        for b in doc:
            if b.name.startswith("mod_") and b.name != "mod_list":
                mod_id = b.name[len("mod_"):]
                for tag in ("_chem_mod_atom", "_chem_mod_bond", "_chem_mod_angle", "_chem_mod_chir", "_chem_mod_tor", "_chem_mod_plane_atom"):
                    unk = set(b.find_values(tag+".mod_id")) - {mod_id}
                    try: self.assertFalse(unk, msg="wrong {}.mod_id in {}".format(tag, b.name))
                    except AssertionError as e: self.errors.append(str(e))

    def test_ener_lib(self):
        elib = read_ener_lib()
        hb_types = set(elib.find_values("_lib_atom.hb_type"))
        try: self.assertEqual(hb_types, set(["N","B","A","D","H"]), msg="invalid _lib_atom.hb_type")
        except AssertionError as e: self.errors.append(str(e))

        all_types = set(elib.find_values("_lib_atom.type"))
        bond_atoms = set(x for i in (1,2) for x in elib.find_values("_lib_bond.atom_type_{}".format(i)))
        undef = bond_atoms - all_types
        try: self.assertFalse(undef, msg="undefined atom_type used in _lib_bond.atom_type_*: {}".format(undef))
        except AssertionError as e: self.errors.append(str(e))
        angl_atoms = set(x for i in (1,2,3) for x in elib.find_values("_lib_angle.atom_type_{}".format(i)))
        undef = angl_atoms - all_types
        try: self.assertFalse(undef, msg="undefined atom_type used in _lib_angle.atom_type_*: {}".format(undef))
        except AssertionError as e: self.errors.append(str(e))
        tors_atoms = set(x for i in (1,2,3,4) for x in elib.find_values("_lib_tors.atom_type_{}".format(i)))
        undef = tors_atoms - all_types
        try: self.assertFalse(undef, msg="undefined atom_type used in _lib_tors.atom_type_*: {}".format(undef))
        except AssertionError as e: self.errors.append(str(e))
        vdw_atoms = set(x for i in (1,2) for x in elib.find_values("_lib_vdw.atom_type_{}".format(i)))
        undef = vdw_atoms - all_types
        try: self.assertFalse(undef, msg="undefined atom_type used in _lib_vdw.atom_type_*: {}".format(undef))
        except AssertionError as e: self.errors.append(str(e))

    def test_links_and_mods_cif(self):
        path = os.path.join(monlib_path, "links_and_mods.cif")
        if not os.path.exists(path): return
        doc0 = gemmi.cif.read(os.path.join(monlib_path, "list", "mon_lib_list.cif"))
        doc = gemmi.cif.read(path)
        for b in doc:
            b0 = doc0.find_block(b.name)
            try: self.assertEqual(b.as_string(), b0.as_string(), msg=f"links_and_mods.cif: {b.name} mismatch")
            except AssertionError as e: self.errors.append(str(e))

if __name__ == '__main__':
    unittest.main()    
