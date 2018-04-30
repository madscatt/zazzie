import sassie.build.pdbscan.pdbscan as pdbscan

pdbfile = 's.pdb'

mol = pdbscan.SasMolScan()
mol.read_pdb(pdbfile)
mol.run_scan()

print 'chain info'
for k,v in mol.chain_info.sequence.iteritems():
    chain_sequence =  v
    print 'chain : ' + k + '\n\tsequence = %s\n' % (chain_sequence,)

print 'segname info'
for k,v in mol.segname_info.sequence.iteritems():
    segname_sequence =  v
    print 'segname : ' + k + '\n\tsequence = %s\n' % (segname_sequence,)
