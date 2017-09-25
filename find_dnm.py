''' Given a list of pedigrees, finds all de novo mutations called by different packages '''


import sys


# where to find results files
result_root = {'penncnv': '/data/NCR_SBRB/big_fake_simplex/penncnv/%s.triocnv',
	       'penncnv_adjusted': '/data/NCR_SBRB/big_fake_simplex/penncnv/%s.adjusted.triocnv'}

ped_list = sys.argv[1]
fid = open(ped_list, 'r')
ped_files = [l.rstrip() for l in fid]
fid.close

ped_names = [p.replace('.ped','').split('/')[-1] for p in ped_files]

for pack, root in result_root.iteritems():
	print 'Looking for De Novo Mutations from %s' % pack
	for f in range(len(ped_files)):
		[kid, mom, dad] = id_family(ped_files[f])
		if pack.find('penn') == 0:
			res = filter_penncnv(root, ped_names[f], kid, mom, dad)


print ped_files
print ped_names
