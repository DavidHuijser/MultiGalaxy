#!/usr/bin/python
import dnest4.classic as dn4
dn4.postprocess(single_precision=True)
import display
import numpy as np


posterior_sample = dn4.my_loadtxt("posterior_sample.txt")
indices = dn4.load_column_names("posterior_sample.txt")


# Extract column by name
#print("Indices")
#print(indices)

#plt.show()
print(indices)
print(type(indices))

print("Length of the dictionanry : ", len(indices))


I = indices.get('indices') 

print("Items of the dictionanry : ", indices.items() )


Z = indices.get('colnames')
print("Colnames: ", Z)
print("Type of Colnames: ", type(Z))
print("Type of items : ", type(I))
print("Type of Posterior  ", type(posterior_sample)  )


index = 'log_a[2]'
print( "Index: ", I['log_a[2]'])
index = I['dim_components']
print( "Index: ", I['dim_components'])
print(posterior_sample[:,int(index)])
print("Print #dim_components:",  posterior_sample[:, I['dim_components']])



#print "dict['Name']: ", dict['Name']
#print "dict['Age']: ", dict['Age']

#plt.hist(posterior_sample[:, indices["num_components"]], 100)

# Extract column by name
#plt.hist(posterior_sample[:, indices[], 100)
#plt.show()
