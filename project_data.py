import numpy

def project_data(original_data, eigenvectors_file):
  print "Projecting data to eigenvector space..."
  
  #loads the eigenvectors
  eigenvectors = numpy.loadtxt(eigenvectors_file, delimiter=",").astype(float)
  
  #formating data
  nentries = len(original_data)
  dimensions = len(original_data[0])
  #gets the sum of entries in each column
  vals_sum = numpy.sum(original_data, axis=0)
  transposed_data = [[0 for i in range(nentries)] for j in range(dimensions)]
  for i in range(dimensions):
    for j in range(nentries):
      transposed_data[i][j] = original_data[j][i]
      
  #makes the projection
  projected_data = [[0 for j in range(len(eigenvectors))] for j in range(nentries)]
  for irow in range(len(eigenvectors)):
    for ientry in range(nentries):
      entry = 0
      for icol in range(len(eigenvectors[0])):
	#FinalData = RowFeatureVector x RowDataAdjusted					#data mean
	entry += eigenvectors[irow][icol]*(transposed_data[icol][ientry] - vals_sum[icol]/float(nentries-1))
      projected_data[ientry][irow] = entry
  
  return projected_data