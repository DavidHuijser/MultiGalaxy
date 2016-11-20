#!/usr/bin/python
from pylab import *

import dnest4.classic as dn4
dn4.postprocess(single_precision=True)
#import display
import numpy as np
import pyfits

import os

import matplotlib.pyplot as plt
#import pylab as P

# load data 
hdulist = pyfits.open('Data/fitim.fits')
data = hdulist[0].data # assuming the first extension is a table
data_shape = data.shape
xsize =  data_shape[0]
ysize =  data_shape[1]

# load sample / posterior 
posterior_sample = dn4.my_loadtxt("sample.txt")
indices = dn4.load_column_names("sample.txt")

# obtain the indices / column names from the object
I = indices.get('indices')                              # I is a dictionary
Z = indices.get('colnames')                             # Z is a list


S =  posterior_sample.shape[0]  #  number of samples
L =  posterior_sample.shape[1]  #  total length of array parameters + hyperparameters + image 

print(" \n ")
print("the dimensions of the image are: ", xsize , " x " , ysize)
print ("Posterior shape: ", posterior_sample.shape)
print ("Length: ", posterior_sample.shape[1])
print ("Sample Size: ", posterior_sample.shape[0])
print ("Starting Position: ", posterior_sample.shape[1] - xsize*ysize)
starting = posterior_sample.shape[1] - xsize*ysize
print(" \n ")

#obtain the maximum number of objects
index = I[Z[2] ]  
temp = posterior_sample[0:,int(index) ]
NOO = temp[0]
print("Max number of objects", NOO)

# obtain starter indices for each type of parameter
# titels[8] = "mag-bar"; titels[9] = "Rout";  titels[10] = "a"; titels[11] = "b"; titels[12] = "q-bar"; titels[13] = "theta-bar";  titels[14] = "box-bar"; 

# print all the names and the first pixel of the imaghe
print(Z[0:52]) 

all_names =  ['sigma', 'dim_components', 'max_num_components', 'hyper_location', 'hyper_scale', 'num_components', 'x', 'x', 'x', 'y', 'y', 'y', 'mag', 'mag', 'mag', 'Re', 'Re', 'Re', 'n', 'n', 'n', 'q', 'q', 'q', 'theta', 'theta', 'theta', 'boxi', 'boxi', 'boxi', 'mag-bar', 'mag-bar', 'mag-bar', 'Rout', 'Rout', 'Rout', 'a', 'a', 'a', 'b', 'b', 'b', 'q-bar', 'q-bar', 'q-bar', 'theta-bar', 'theta-bar', 'theta-bar', 'box-bar', 'box-bar', 'box-bar']

print(size(all_names), len(all_names) ) 
print(posterior_sample[:,6:9])

all_index = np.empty(len(all_names))

num_comp_index = Z.index('num_components')
total_objects = np.sum(posterior_sample[:,num_comp_index])


print("total number of objects ", total_objects)
print("Objects ", posterior_sample[:,num_comp_index])


#for i in range(6, len(all_names)):
for i in range(0,len(all_names)):
      print(all_names[i],  " ", Z.index(all_names[i]))
      this_index = Z.index(all_names[i])

      if i > 5:
          subset = posterior_sample[:,this_index:this_index+NOO]                   # shape (samplesize, max_num_components)
          index = [subset != 0]
      else:
          subset = posterior_sample[:,i]                                           # shape (samplesize, max_num_components)
          index = [subset != 0]
      print(subset.shape)               

      # in case of magnitude rescale to prior-distribution
      if ((Z[i] == 'mag') or (Z[i]== 'mag-bar')):
         print("Rescale Magnitude")
         for k in range(0, S):
            location_index = posterior_sample[k,3]              # 3 
            scale_index =  posterior_sample[k,4]                # 4 
            for j in range(0,int(NOO)):
                  subset[k,j] = (subset[k,j]-location_index )/scale_index



#         print(location_index.shape, scale_index.shape, subset.shape)

#         new_subset = help_subset[subset != 0]           
#      else:
#         new_subset = subset[subset != 0]           

     # params[i] = 5*rt(2) + 25;   // shape=2, location=25, scale=5

      new_subset = subset[index]   


      title = 'Histogram ' +all_names[i] 
      plt.title(title) 
      plt.xlabel(all_names[i])

      plt.hist(new_subset, alpha=0.5, bins=20)
      filename = 'test_prior_histogram_' +all_names[i] +'.png'
      plt.savefig(filename)      
#      if (i == 0):
#           i = 2

#      if (i > 5):
#           i = i+ NOO

#      new_subset = subset[subset != 0]
#      print((new_subset.shape)
      # make histogram for each 
#      P.figure()
 #     n, bins, patches = P.hist(subset[subset != 0], 50, normed=1, histtype='stepfilled')
  #    P.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
   #   P.show()
	


#      sub = []
#      this_index = Z.index(all_names[i])
#      for j in range(0, S):
#         gap = posterior_sample[j,Z.index(all_names[5])] 
#         print(gap)
#         sub.append(posterior_sample[j,this_index:this_index+gap])
#      print(sub)
      
      #make a subset of all xvalues using a list
      


#all_index = Z.index(all_names[0])
#y_index = Z.index('y')

#mag_index = Z.index('mag')
#Re_index = Z.index('Re')

#n_index = Z.index('n')
#q_index = Z.index('q')

#theta_index = Z.index('theta')
#boxi_index = Z.index('boxi')

#mag_bar_index = Z.index('mag-bar')
#Rout_index = Z.index('Rout')

#print(a_index = Z.index('a'))
#print(b_index = Z.index('b'))

#print(q_bar_index = Z.index('q-bar'))
#print(theta_index = Z.index('theta-bar'))

#print(box_bar_index = Z.index('box-bar'))
#print(sigma_index = Z.index('sigma'))


#print(hyper_scale_index = Z.index('hyper-scale'))
#print(hpyer_location_index = Z.index('hpyer-location_'))


#print(x_index, y_index)



## posterior_sample[:, I['']])
#print("Element 1", Z[0]," ")  # sigma
#print("Element 2", Z[1]," ")  # dim comp
#print("Element 3", Z[2]," ")  # max num comp

#print("Element 4", Z[3]," ")  # hyper location
#print("Element 5", Z[4]," ")  # hyper scale
#print("Element 6", Z[5]," ")  # num comp

#print("Element 7", Z[6]," ")  # x
#print("Element 8", Z[7]," ")  # x
#print("Element 9", Z[8]," ")  # x


#print("Element 10", Z[9]," ")  # x
#print("Element 11", Z[10]," ")  # x
#print("Element 12", Z[11]," ")  # x

#print(" \n ")

#print("print keys", indices.keys())

#print("print keys", I.keys())
#print("print keys", type(Z))
#print("print keys", type(I))  #

##print("print values", indices.values() )
#print( indices.get("x", "none"))
#print( Z[5])
#print( I.get(Z[5], "none"))



#print(" \n ")
#print("Element 1", Z[0]," ")  # sigma
#print("All the sigmas",  posterior_sample[:,0 ])           # non-flexible
#print("All the sigmas",  posterior_sample[:,I==Z[0] ])     # flexible
#print("All the sigmas",  posterior_sample[:,I=='sigma'])   # even more flexible

#print("Element 2", Z[1]," ")  # sigma
#print("All the dim_components",  posterior_sample[:,1 ])           # non-flexible
#print("All the dim_components",  posterior_sample[:,I[Z[1]] ])     # flexible
#print("All the dim_components",  posterior_sample[:,I['dim_components']])   # even more flexible


#print("Element 6", Z[5]," ")  # sigma
#print("All the dim_components",  posterior_sample[:,5 ])           # non-flexible
#print("All the dim_components",  posterior_sample[:,I[Z[5]] ])     # flexible
#print("All the dim_components",  posterior_sample[:,I['num_components']])   # even more flexible


#print("Element 7-9", Z[6]," ")  # sigma
#print("All the dim_components",  posterior_sample[:,6:9 ])           # non-flexible
##print("All the dim_components",  posterior_sample[:,I[Z[6:9]] ])     # flexible
##print("All the dim_components",  posterior_sample[:,I['x']])   # even more flexible



##print("Element 6", Z[5]," ")  # sigma
##print("All the num comp",  posterior_sample[:,6 ])           # non-flexible
##print("All the num comp",  posterior_sample[:,I==Z[5] ])     # flexible
##print("All the num comp",  posterior_sample[:,I=='num_components'])   # even more flexible



##print("Element ", Z[6]," ")  # x-vales
##print("All the x",  posterior_sample[:, ])           # non-flexible
##print("All the x",  posterior_sample[:,I[Z[6]] ])     # flexible
##print("All the x",  posterior_sample[:,I=='x'])   # even more flexible


##index = np.where(Z=='sigma')
##print("index: ", index," ")
##print("index: ", index," ",Z[index]," ", I[index] )




##starter_point = 6
##subset = posterior_sample[:,I[Z[7]] ]
##print("Subsetshape", subset)
##print("Element 7", Z[6:6+3]," ")  # x

##print("X coordinates",posterior_sample[:,0:51] ,"\n")




##index = I['num_components']
##print( "Index: ", I['num_components'])
##print(posterior_sample[:,int(index)])



#index = I[Z[6] ]
#temp = posterior_sample[0:,int(index) ]
#NOOx = temp[0]
#print("X", NOOx)
#print("X", NOOx+NOO)
##print(posterior_sample[:,int(index)+1:int(index)+1+NOO])


###img = posterior_sample[9, starting:starting+xsize*ysize].reshape((xsize, ysize))            # Need to edit this
####img = posterior_sample[0, starting:starting+xsize*ysize].reshape((ysize, xsize))            # Need to edit this  
###print("Shape of image is", type(img))

###subplot(1, 3, 2)

###fig = plt.figure(figsize=(8,6))
###ax = fig.add_subplot(111)
###fig.subplots_adjust(top=0.85)
###imshow(img, cmap='copper', interpolation='nearest')
###filename = 'first_data.png'
###plt.savefig(filename)
###filename = 'first_data.pdf'
###plt.savefig(filename)
###filename = 'first_data.eps'
###plt.savefig(filename)


### Extract column by name
###print("Indices")
###print(indices)

###plt.show()
###print(indices)
###print(type(indices))

###print("Length of the dictionanry : ", len(indices))


##I = indices.get('indices') 

###print("Items of the dictionanry : ", indices.items() )


##Z = indices.get('colnames')
###print("Colnames: ", Z)
##print("Type of Colnames: ", type(Z))
##print("Type of items : ", type(I))
##print("Type of Posterior  ", type(posterior_sample)  )


###index = 'log_a[2]'
###print( "Index: ", I['log_a[2]'])
##index = I['num_components']
##print( "Index: ", I['num_components'])
##print(posterior_sample[:,int(index)])



###print("Print #num_components:",  posterior_sample[:, I['num_components']])
###print("Print #log_box-bar[2]:",  posterior_sample[:, I['log_box-bar[2]']])
###<type 'dict'>
###('Colnames: ', ['sigma', 'dim_components', 'max_num_components', 'hpyer_location_', 'hyper_scale', 'num_components', 'x', 'x', 'x', 'y', 'y', 'y', 'mag', 'mag', 'mag', 'Re', 'Re', 'Re', 'n', 'n', 'n', 'q', 'q', 'q', 'theta', 'theta', 'theta', 'boxi', 'boxi', 'boxi', 'mag-bar', 'mag-bar', 'mag-bar', 'Rout', 'Rout', 'Rout', 'a', 'a', 'a', 'b', 'b', 'b', 'q-bar', 'q-bar', 'q-bar', 'theta-bar', 'theta-bar', 'theta-bar', 'box-bar', 'box-bar', 'box-bar'

###print("Print #log_box-bar[2]:",  posterior_sample[:, I['']])



