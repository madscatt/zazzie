import sasmol.system as system



pdb_files = ["1E3B_label.pdb","2ZIL_lysozyme_Au_label.pdb", "DNA10Au4Nan.pdb", "DNA10AuTbNan.pdb", "1E3B_label.pdb", "2ZIL_lysozyme_Au_label.pdb", "DNA10Au4Nan.pdb", "DNA10AuTbNan.pdb"]

for this_pdb in pdb_files:

    m = system.Molecule(this_pdb)

    print("%s has %i atoms\n" % (this_pdb, m.natoms()))

    
