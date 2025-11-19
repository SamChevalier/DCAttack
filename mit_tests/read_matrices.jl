using HDF5

# read the matrices directly
testdata_file = "./mit_tests/ABc_data.h5"
fid           = h5open(testdata_file, "r")
A             = read(fid, "A")  
B             = read(fid, "B") 
c             = read(fid, "c")
close(fid)
