from PIL import Image
import numpy as np

#user define x-range of image, how far it is from lens (in z direction)
#and how fast it's moving in z-direction
XMin=-10.0
XMax=10.0
YOffset=0.0 # how far to offset Y, 0 means centered
InitZ=-100.0
InitVz=1

#load image and convert to array
img = Image.open('ladybug.jpg')
ar = np.array(img)

#Use # of pixels to define y-range of image, keeping aspect ratio
StepSize=ar.shape[1]/(XMax-XMin)
YMax=ar.shape[0]/(2*StepSize)+YOffset
YMin=-ar.shape[0]/(2*StepSize)+YOffset

#define x,y,z positions and velocitites for each pixel 
xpos=np.linspace(XMin,XMax,num=ar.shape[1])
xvel=np.zeros(len(xpos))
ypos=np.linspace(YMin,YMax,num=ar.shape[0])
yvel=np.zeros(len(ypos))
zpos=InitZ
zvel=InitVz

x,y,z=np.meshgrid(xpos,ypos,zpos,indexing='ij')
vx,vy,vz=np.meshgrid(xvel,yvel,zvel,indexing='ij')
R=ar[:,:,0]
G=ar[:,:,1]
B=ar[:,:,2]
