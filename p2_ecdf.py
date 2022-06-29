import pickle

probe_rpms = pickle.load(open('probe_rpms.pickle','rb'))
probe_ecdfs = pickle.load(open('probe_ecdfs.pickle','rb'))

print(probe_rpms)
print(probe_ecdfs)
