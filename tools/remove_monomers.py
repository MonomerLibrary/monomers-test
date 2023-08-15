import gemmi
if __name__ == "__main__":
    import sys
    lst_in = sys.argv[1]
    names = set(sys.argv[2:])
    
    doc = gemmi.cif.read(lst_in)
    block = doc.find_block("comp_list")
    table = block.find("_chem_comp.", ["id"])
    todel = []
    todel_files = []
    for i, r in enumerate(table):
        mon = r.str(0)
        if mon in names:
            todel.append(i)
            todel_files.append("{}/{}.cif".format(mon[0].lower(), mon))
    for i in reversed(sorted(todel)):
        table.remove_row(i)

    print("git rm {}".format(" ".join(todel_files)))
    doc.write_file(lst_in)
