import tensorflow as tf


a=tf.constant([1,3,5,7,8,5,4,8,9,5])
b=tf.constant([1,3,5,7,8,5,4,8,9,5])
result=tf.multiply(a,b)

config=tf.ConfigProto(log_device_placement=True)
#print(config)
with tf.Session(config=tf.ConfigProto(log_device_placement=True)) as sess:
  output = sess.run(result)
  print(output)
sess.close()