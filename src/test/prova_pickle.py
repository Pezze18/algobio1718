import pickle

fileObject = open("b", 'wb')
pickle.dump(fileObject, 10)
pickle.dump(fileObject, 2)
fileObject.close()