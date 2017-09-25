''' Given a list of pedigrees, finds all de novo mutations called by different packages '''


import sys


# where to find results files
result_root = {'penncnv': '/data/NCR_SBRB/big_fake_simplex/penncnv/'}

ped_list = sys.argv[1]
fid = open(ped_list, 'r')
ped_files = [l.rstrip() for l in fid]
fid.close

ped_names = [p.replace('.ped','').split('/')[-1] for p in ped_files]

for pack, root in result_root.iteritems():
	print 'Looking for De Novo Mutations from %s' % pack

print ped_files
print ped_names
