# -------------------------------------- INFO --------------------------------------
'''
Date: 2023/10/23
Version: 1.0
Main Author: Lorenzo Boccia lorenzo.boccia@unina.it - University of Naples Federico II
Co-author: Gabriele Delogu gabriele.delogu@unitus.it - Tuscia Univesity
Co-author: Miriam Perretta miriam.perretta@unina.it - University of Naples Federico II
Website: www.larp.unina.it
Use is free of charge

'''
# -------------------------------------- AIM & NOTES--------------------------------------

'''
Prisma Tool is a Python script. It can be used to open HDF5 data from the PRISMA
satellite mission and convert them to Tiff. If necessary, the software can also be used
to modify the coordinates of the satellite image vertices.
It seems to work for L1 and L2D data from PRISMA (not tried for L2B or L2C).
The program opens the HDF5 file, saves the metadata and the tree of the HDF5 file,
and displays some values on the screen and generates the a Tiff image for each band.
Before starting this program, please read the Word document -Prisma_Tool_installation_guide.doc-.
should be consulted. It contains the installation guide.
If you develop a new version, please send it to lorenzo.boccia@unina.it
The code was developed in Italian.
Therefore the names of some functions or comments are still in Italian.
We will upload a new version fully translated into English as soon as possible!
'''
  
# -------------------------------------- IMPORT LIBRARIES --------------------------------------

#  h5py is the library for managing files in the HDF5 Hierarchical Data Format
import  h5py

# The os library is used to open directories and handle files.
# For example to know the directory or to create it
import os

# from the sys library imports the function to terminate execution
import sys
from sys import exit

# enable tkinter (activates the GUI tkinter 
# is the library that includes the window manager filedialog,
# to open the window in which to choose files to open)
from tkinter import filedialog

# now we import the whole tkinter library with *
# from now on commands (functions and classes) should not be preceded by tkinter.
# For example I don't have to write tkinter.title but simply title.
# having imported the entire library perhaps no longer needing to import filedialog
from tkinter import *

# We import under the name np the Numpy library
# it provides many mathematical functions, including matrices
# and is used to run raster but also to handle matrices.
import numpy as np

# import maths to have sine and cosine arcotangent
import math

# Import rasterio which is the 'raster input output' library that actually runs
# the GDAL library which is in C++ we need to export images in geotiff format. 
# Rasterio needs to be imported from the terminal;
# command: python -m pip install rasterio as in the Word instructions. 
import rasterio.transform
import affine

# the pyproj library is a coordinate converter (based on the EPSG code) from which we import the Transformer function
from pyproj import Transformer

# the matplotlib pyplot called plt is the library for displaying rasters on screen in python
import matplotlib.pyplot as plt

# -------------------------------------- MATRIX DIMENSIONING --------------------------------------

# Let's dimension some matrices, defining them as numpy objects so that we can work much faster.
# L1 images are smaller, but we need to dimension them knowing that L2D images are not 1000*1000 or 6000*6000, but larger! 
# In L2D data, the frame is large and contains the image, which is rotated,
# so we don't have to dimension 1000*1000 but 2000*2000.
# Maybe 1500*1500 would be enough for the diagonal of the square.
# For the Pancromatic, 6000*6000 is not enough, but 9000*9000 would be enough.

vnir =np.empty([2000,65,2000])
swir =np.empty([2000,172,2000])
imm_banda = np.empty([2000,2000])
pancromatico=np.empty([12000,12000])
latitudine_pan=np.empty([12000,12000])
longitudine_pan=np.empty([12000,12000])

# Vnir and swir are cubic sextinate. Pancromatic to the hyperspectral image.
# Imm_band is allocated to the single hyperspectral image to be saved in the geotiff.
# Latitude and Longitude are allocated to the pixel coordinates
# L1 cube 1000*1000 Vnir and Swir, 6000*6000 panchromatic (page 98-261 prisma user manual)

# -------------------------------------- FUNCTIONS --------------------------------------

def estrai_banda(matrice, banda_numero):
    '''matrix can be swir or vnir, the size is 1000*1000 if it is an L1. band_number can be 0-65 for VNIR 0-172 for SWIR'''
    estratto = np.empty([2000, 2000])
    estratto = matrice[:, banda_numero, :]
    # in python indexes e.g. 2:18 from second to eighteenth.
    # With ":" we mean all. With "2 :" we mean from the second to all
    return estratto

def proiezione (lat, long):
    '''Enter with lat and long matrices and exit with the epsg projection number.
    Basically we enter with latitude and longitude and try to figure out which UTM projection
    we should use based on the hemisphere and the zone'''
    # The epsg 4326 is wgs84
    # if we are west of 12° we must use UTM 32 (epsg 32632) otherwise UTM 33 (epsg 32633)
    lat_min=np.min(lat)
    long_min=np.min(long)
    terzacifra= (32600)
    # northern hemisphere begin with 326xx, southern begin with 327xx
    if lat_min < 0:  
        terzacifra= (32700)
    proiezione = terzacifra + 30 + int((long_min/6+1))
    return proiezione

def vertici (lat, long, proiezione, iper,tipo_prodotto,sequenza,ri_georeferenzia,deltageoref):
    
    '''as input the matrix of latitudes, longitudes and the projection
    as output a table with 6 rows (i,J, lat, long, xcoord and ycoord),
    4 columns of the 4 vertices and the two dimensions in m of the GSD.
    ri_georeferenzia is "1" if re_georeferenced is needed.
    sequence are the headers for example "top- right, northeast, ...".
    deltageoref are the vertices of the panchromatic which are stored to calculate those of the hyperspectral.'''
    tabella=np.empty([6,4])
    '''Since the first two rows of the table are indexes, newer versions of numpy use only integers as indexes
    (or at least that's what happened to us). Python actually allows mixed lists of integers and reals, but to
    gloss over the problem we created an array of integers ('spigoli') with indices and then added it to the array
    in the first two rows.'''
    spigoli=np.empty([6,4], dtype=int)
    # lat is the latitude of the n*m pixels where: n is across (East) track and m is along (North),
    # pixels are 6000 indexed from 0 to 5999
    # in the L2D case, the first line (0.0 to 7283.0) goes from the Northwest edge to the Southwest edge and is
    # therefore in the direction of the satellite (along track) going from North to South 
    pixel_along,pixel_across =lat.shape
    # we assign to "spigoli[,]" because it is declared of integers, then we copy them into "tabella[,]"
    spigoli[0,0]=0; spigoli[0,1]=pixel_along-1; spigoli[0,2]=0; spigoli[0,3]=pixel_along-1
    spigoli[1,0]=0; spigoli[1,1]=0; spigoli[1,2]=pixel_across-1; spigoli[1,3]=pixel_across-1
    
    # In the case of L2D these are rows and columns of the vertices of the outer rectangle, not the image scene!
    # To convert from geographical to UTM there are two procedures in which the order of lat and log is reversed.
    # Procedure 1 which works but gives warnings:
    # lonlat=pyproj.Proj(init='epsg:4326')
    # utm=pyproj.Proj(init='epsg:32633')
    # est,north=pyproj.transform(lonlat,utm,longitudine_pan[0,0],latitudine_pan[0,0])
    # procedure 2 which works better but is inverted lat and long
    # est1,north1=pyproj.transform(4326,32633,latitudine_pan[0,0],longitudine_pan[0,0])
    for cont_j in range (0,4):
        tabella[2, cont_j] = lat[spigoli[0,cont_j],spigoli[1,cont_j]]
        tabella[3, cont_j] = long[spigoli[0, cont_j], spigoli[1, cont_j]]  
        if ri_georeferenzia== 1 and iper==0:
            if tipo_prodotto=="PRS_L2D_STD":
                print ("It's an L2D product. To Georeferencing it you have to refer to the external frame")
                print ("you'll need the new coordinates of the center of the pixel for each corner of the pan FRAME")
            else:
                print ("It's an ", tipo_prodotto," product. To Georeferencing it you have to refer to immage")
                print ("you'll need the new coordinates of the center of the pixels for each corner of the pan IMMAGE")
                print ("IMAGE corresponds to vertices of square with non-null values, not vertices of outer frame.")
            # The values in deltageoref are useful to avoid having to retype data a second time
            deltageoref[2,cont_j]=tabella[2, cont_j]
            deltageoref[3,cont_j]=tabella[3, cont_j]
            print ('Corner ', sequenza[cont_j],' of the image ')
            print ('original LATitude  was                          ' , tabella[2,cont_j])
            print ('new LATitude of corner is ',sequenza[cont_j], end="\t")
            tabella[2,cont_j]= float(input (' '))
            print ('Original LONGitude was                          ' , tabella[3,cont_j])
            print ('New LONGitude of the corner is ',sequenza[cont_j] , end="\t")
            tabella[3,cont_j]= float(input (' '))
            deltageoref[2,cont_j]=deltageoref[2,cont_j]-tabella[2, cont_j]
            deltageoref[3,cont_j]=deltageoref[3,cont_j]-tabella[3, cont_j]
            # this way we have the deltas of georeferencing preserved
        elif ri_georeferenzia== 1 and iper==1:
            # means we are in the hyperspectral; we maintain the differences in lat and long
            tabella[2, cont_j] = tabella[2, cont_j]-deltageoref[2, cont_j]
            tabella[3, cont_j] = tabella[3, cont_j]-deltageoref[3, cont_j] 
        # with 'cambio' we create the function to transform from one system to another via the Pyproj library
        cambio=Transformer.from_crs(4326, proiezione)
        tabella[4, cont_j], tabella[5, cont_j]= cambio.transform(tabella[2, cont_j], tabella[3, cont_j])
    # at this stage we join the 'spigoli' table of integers to the 'tabella' in the first two rows so that the table contains
    # in each column the indices, the lat and long of the vertices and the UTM North and East of the vertices
    for cont_j in range (0,4):
        tabella[0,cont_j] = spigoli[0,cont_j]
        tabella[1,cont_j] = spigoli[1,cont_j]
    # Dim along is the pixel size of the L2D image in the NO SO direction
    # Dim across is the pixel size of the L2D image in the NO NE direction
    dim_across = ((((tabella[5,1]-tabella[5,3]))**2+((tabella[4,1]-tabella[4,3]))**2)**0.5)/(pixel_across-1)
    dim_along = ((((tabella[4,0]-tabella[4,1]))**2+((tabella[5,0]-tabella[5,1]))**2)**0.5)/(pixel_along-1)
    # (pixel_along - 1) because if there are n points there are n-1 intervals (e.g. 7283 if pixel_along is 7284).
    # in pixels along has the shape of the latitude matrix and thus the number of pixels, not the greatest index
    return tabella,dim_across, dim_along, deltageoref

def azimut (punto1,punto2):
    # as input there are the UTM x,y coordinates of two vertices calculate azimuth1-2)
    if abs(punto2[1] - punto1[1]) > 0.001:
        azimut_radian=math.atan ((punto2[0]-punto1[0])/(punto2[1]-punto1[1]))
    else:
        azimut_radian = math.pi / 2
    azimutlato=math.degrees(azimut_radian)+180
    return azimutlato

def stampa_su_file(cosa_da_stampare):
    # This seemingly useless function is required to be able to use the h5py visit method,
    # which must have a function as an argument.
    file_da_scrivere.write(cosa_da_stampare)
    file_da_scrivere.write('\n')

def percorsi_nomi(nome_completo_file):
    # identifying the directory used
    # there is probably an easier way
    last_slash = 0
    numero_carattere = len(nome_completo_file)
    while last_slash <= 0:
        if nome_completo_file[numero_carattere - 1:numero_carattere] == "/":
            last_slash = numero_carattere - 1
        numero_carattere = numero_carattere - 1
    last_slash = last_slash + 1
    posizione_input = nome_completo_file[0:last_slash]
    nome_input = nome_completo_file[last_slash:]
    # last_slash is the location of the last / in the file name
    # position_input is the directory of the file on disk
    # input_name is the name of the source
    # filenames L1 and L2 have different formats (Prisma user guide p.109 vs. 234).
    # I need the position of the fourth underscore which in each case is followed by the year the month and the day etc.
    fourth_ = 0
    vero = 0
    numero_carattere = len(nome_input)
    while fourth_ <= 0:
        if nome_input[numero_carattere - 1:numero_carattere] == "_":
            vero = vero + 1
            if vero == 2:
                fourth_ = numero_carattere
        numero_carattere = numero_carattere - 1
    data = nome_input[fourth_:fourth_ + 4] + "_" + nome_input[fourth_ + 4:fourth_ + 6] + "_" + nome_input[fourth_ + 6:fourth_ + 8]
    # print("data", data)
    
    codice=nome_input[3:7]
    # print (codice[3:])
    
    if codice[3:] != '_':
        # is L2D product
        tipo_prod = codice[1:4]
        codice=codice+'_'
    else:
        tipo_prod = codice[1:3]
    posizione_uscita = posizione_input + "GeoTIFF" +codice+ data + "/"
    print('Recognized', tipo_prod)
    return posizione_input, nome_input, posizione_uscita, tipo_prod, codice, data


def affinity (tipo_prodotto,tab_coord,x_pixel,y_pixel):
    if tipo_prodotto=="PRS_L1_STD" or tipo_prodotto== "PRS_L2B_STD" or tipo_prodotto== "PRS_L2C_STD":
        ne= [tab_coord[4,0],tab_coord[5,0]]
        no= [tab_coord[4,1],tab_coord[5,1]]
        se= [tab_coord[4,2],tab_coord[5,2]]
        so= [tab_coord[4,3],tab_coord[5,3]]
        base=tab_coord[0,1]*x_pixel
        altezza=tab_coord[0,3]*y_pixel
    else:
        no= [tab_coord[4,0],tab_coord[5,0]]
        so= [tab_coord[4,1],tab_coord[5,1]]
        ne= [tab_coord[4,2],tab_coord[5,2]]
        se= [tab_coord[4,3],tab_coord[5,3]]
        base=tab_coord[1,2]*x_pixel
        altezza=tab_coord[0,3]*y_pixel
    # to find out how much we need to rotate we consider the diagonal of the image and the diagonal ne_so
    azimut_ne_so=azimut(ne,so)
    azimut_ne_no=azimut(ne,no)+180
    azimut_ne_se=azimut (ne,se)
    so_immage=[ne[0]+base,ne[1]+altezza]
    azimut_immage=azimut(ne,so_immage)
    # print ("azimut_ne_so",azimut_ne_so)
    # print ("azimut_immage",azimut_immage)
    alfa= azimut_ne_so - azimut_immage
    # print ("alfa", alfa)
    # one must calculate the shear of the two side shear1 and shear2
    shear1=270-azimut_immage-(azimut_ne_no-azimut_ne_so)
    shear2=(azimut_immage-180)-(azimut_ne_so-azimut_ne_se)
    #print((azimut_immage-180),(azimut_ne_so-azimut_ne_se))
    #print ("shear1,shear2", shear1,shear2)
    #print ("azimut_ne_no",azimut_ne_no)
    #print ("azimut_ne_se",azimut_ne_se)
    alfa=alfa+90
    #print ("alfa", alfa)
    #input("controllo alfa")
    if tipo_prodotto=="PRS_L1_STD" or tipo_prodotto== "PRS_L2B_STD" or tipo_prodotto== "PRS_L2C_STD":
        #print ("alfa", alfa)
        #print("alfa+270", alfa+180)
        rotazione = affine.Affine.rotation(alfa + 180)
        a, b, c, d, e, f, g, h, i = rotazione
        #print(a, b, c, d, e, f, g, h, i, "rot a,b,c,d,e,f,g,h,i")
        # ROTAZIONE è riferita al lato NE NO .
        scala = affine.Affine.scale(x_pixel, y_pixel)
        a, b, c, d, e, f, g, h, i = scala
        #print(a, b, c, d, e, f, g, h, i, "scale a,b,c,d,e,f,g,h,i")
        scala_rot = rotazione * scala
        a, b, c, d, e, f, g, h, i = scala_rot
        #print(a, b, c, d, e, f, g, h, i, "scale rot a,b,c,d,e,f,g,h,i")
        # determination of the pivot:
        # the pivot is in the centre of the pixel and not in the image corner, so it makes sense not
        # to do this calculation and leave the pivot in the centre of the pixel
        # angolo=math.radians(45-(alfa-180))
        # xpivot_pan=tabella_pan[4,0]+x_pan_pixel/(2**0.5)*math.sin(angolo)
        # ypivot_pan=tabella_pan[5,0]+x_pan_pixel/(2**0.5)*math.cos(angolo)
        xpivot=tab_coord[4,0]
        ypivot=tab_coord[5,0]
        traslazione=affine.Affine.translation(xpivot,ypivot)
        a, b, c, d, e, f, g, h, i = traslazione
        scala_rot = rotazione * scala
        a, b, c, d, e, f, g, h, i = scala_rot
        # print(a, b, c, d, e, f, g, h, i, "scala rot a,b,c,d,e,f,g,h,i")
        # print ("dist_shear= ",dist_shear," alfa=", alfa)
        deformazione = affine.Affine.shear(+shear1, +shear2)
        scala_rot_shear = scala_rot * deformazione
        a, b, c, d, e, f, g, h, i = scala_rot_shear
        # a b c we assign them to negative to have reflection
        # put the pivot in c and f
        a = -a
        b = -b
        d = d
        e = e
        # c=tab_coord[4,0]+x_pixel/2
        # f=tab_coord[5,0]+y_pixel/2
        c=tab_coord[4,0]
        f=tab_coord[5,0]
        # print(a, b, c, d, e, f, g, h, i, "new parameters a,b,c,d,e,f,g,h,i")
    else:
        a = x_pixel
        b = 0
        # c=xpivot
        #c = tab_coord[4, 0]-x_pixel/2
        c = tab_coord[4, 0]
        d = 0
        e = -y_pixel
        # f=ypivot
        #f = tab_coord[5, 0]+y_pixel/2
        f = tab_coord[5, 0]
    trasformazione = a,b,c,d,e,f
    return trasformazione

# -------------------------------------- MAIN --------------------------------------
# -------------------------------------- Management of tree and metadata --------------------------------------

print('Author: Lorenzo Boccia lorenzo.boccia@unina.it  - Univesity of Naples "Federico II" - www.larp.unina.it')
print('Coauthor: Gabriele Delogu gabriele.delogu@unitus.it  - Tuscia Univesity - www.larp.unina.it')
print ('The program will generate GeoTiff on the basis of the HDF5 files of PRISMA satellite L1 or L2D product')

# We open the root window of the tkinter library tk(), which must be destroyed at the end.
# Tk() is the function of the tkinter module that opens a small window.
radice = Tk()
# We open the dialogue box (which is logically nested in the root) to ask or search for the file name.
# In the window with the file name, I open a dialogue box and ask for the file name
nome_file =  filedialog.askopenfilename(initialdir = "/",title = "Select file",filetypes = (("HDF5","*.he5"),("all files","*.*")))
# print (nome_file)
radice.destroy()
# When we received the filename we close the radice.destroy
# "root.mainloop() =" waits for the value of the tkinter window, but in this case it is not needed.
# It is useful when we "write on screen", for example.
posizione,nome_file_sorgente, posizione_output, tipo_prodotto,codice, data = percorsi_nomi(nome_file)
# product_type is L1 or L2D, code is _L2D_ or _L1_

# we open the folder with an os library command
# os.mkdir(location + "GeoTIFF/"+data) would open it, but if it already exists it would give an error.
# with os.makedirs( name,exist_ok=false) we can manage the existing directory
os.makedirs(posizione_output, exist_ok=True)
nome_file_metadato= posizione_output + "Metadato"+data+".txt"
file_da_scrivere=open(nome_file_metadato, 'w')
# open the file we call 'immagine' with the h5py library
# open the existing read-only file 'r'. If you need to edit an
# existing file use 'r+' to read write or create use 'a'.
immagine = h5py.File(nome_file, 'r')

# reading name and value for root attributes (metadata contained in HDF5 root)
# in " immagine.attrs" there is a list or tuple indicating how many components each attribute group has
# "image.attrs" is a group similar to a dictionary
# we need the vectors with the wavelengths of the cw bands and the active bands flags
# then I need to get them from the metafile. For now let's create (declare) four vectors and call them list_..._..
list_cw_swir=[]
list_cw_vnir=[]
list_cw_swir_flags=[]
list_cw_vnir_flags=[]
# Now we need the vectors with the wavelengths of the bands and the vectors with the flags of the active bands
for attributo in immagine.attrs:
    # as "attributo" increases, the variable "contenuto" scrolls the imaggine.attrs dictionary (the metadata) and puts strings in "contenuto"
    contenuto = str(immagine.attrs[attributo])
    # Carica i 4 vettori
    if attributo =="List_Cw_Swir":
        list_cw_swir=(immagine.attrs[attributo])
    if attributo =="List_Cw_Swir_Flags":
        list_cw_swir_flags=(immagine.attrs[attributo])
    if attributo =="List_Cw_Vnir":
        list_cw_vnir = (immagine.attrs[attributo])
    if attributo =="List_Cw_Vnir_Flags":
        list_cw_vnir_flags=(immagine.attrs[attributo])
    # We need to know what type of product is: L1 or L2. The information is in the Product_ID field of the HDF5 file.
    # When I download it to product_type I get an array of numpy.bytes_ I convert it to asci with utf-8
    if attributo=="Product_ID":
        tipo_prodotto=(immagine.attrs[attributo])
        tipo_prodotto=tipo_prodotto.decode('UTF-8')
    file_da_scrivere.write(attributo)
    if attributo=="Epsg_Code":
    # L2d images already have the EPSG code in the metadata. L1 L2c images do not have it and so it remains empty.
        codice_epsg=(immagine.attrs[attributo])
        print (codice_epsg)
        
    # Use \n to line up, otherwise it will write everything on one line.
    file_da_scrivere.write('\n')
    file_da_scrivere.write(contenuto)
    file_da_scrivere.write('\n')
file_da_scrivere.close()

# We use the .visit(callable) method on the image to print out information.
# The visit method takes a function as an argument and not a variable
# we define the function "stampaaschermo(cosa_da_stampare)" which we will need later
# the function results in a print command for example:
# def stampa_a_schermo(cosa_da_stampare): print(cosa_da_stampare)
# the .visit method is a method we get from h5py. We need it to print the tree of the image file
# visit takes a function (in our case "print_screen") as an input variable
# and applies it to all objects in the image group. Basically in the visit we tell it to do
# the screen print of all objects in the groups and subgroups
# not the values (the datasets) but the dataset names (the tree)
# immagine.visit(stampa_a_schermo)
# the tree is essential to know where the information we are searching for is located!
nome_file_tree= posizione_output + "Tree"+data+".txt"
file_da_scrivere=open(nome_file_tree, 'w')

# stampa_su_file is a function defined to write to an open file.
# This allows us to use h5py's visit method to print the tree
immagine.visit(stampa_su_file)
file_da_scrivere.close()

# -------------------------------------- Reading the HIS Cube --------------------------------------
s1 = '/HDFEOS/SWATHS/'
s2 = tipo_prodotto[:-3]
s3 = 'HCO'
s4 = 'PCO'
s5 = '/Data Fields/'
s6 = '/Geolocation Fields/'
vnir = immagine[s1 + s2 + s3 + s5 + 'VNIR_Cube']
swir = immagine[s1 + s2 + s3 + s5 + 'SWIR_Cube']
pancromatico = immagine[s1 + s2 + s4 + s5 + 'Cube']
latitudine_pan = immagine[s1 + s2 + s4 + s6 + 'Latitude']
longitudine_pan = immagine[s1 + s2 + s4 + s6 + 'Longitude']
# print(pancromatico.shape,"\n",latitudine_pan.shape,"\n",longitudine_pan.shape)

if tipo_prodotto=="PRS_L1_STD":
    # WARNING: HCO and PCO are coregistered data (radiance) and not digital numbers (PRISMA manual p. 94)
    # Caution For L1 coregistered images the coordinates are those of the VNIR (PRISMA manual p. 93)
    # to read SWIR & VNIR datacubes I see where SWIR_Cube is from the tree and read it in swir
    # for L1 type data the file name is suffixed with _VNIR
    latitudine = immagine[s1+s2+s3+ s6+'Latitude_VNIR']
    longitudine = immagine[s1+s2+s3+ s6+'Longitude_VNIR']
    print('The product type is:', tipo_prodotto)
    
elif tipo_prodotto== "PRS_L2D_STD" or tipo_prodotto== "PRS_L2B_STD" or tipo_prodotto== "PRS_L2C_STD":
    latitudine = immagine[s1 + s2 + s3 + s6 + 'Latitude']
    longitudine = immagine[s1 + s2 + s3 + s6 + 'Longitude']
    print('The product type is: ', tipo_prodotto)
    
else:
    print('Something strange or not yet manageable in',tipo_prodotto)
    exit()
    
# list_cw_vnir=image['/HDFEOS/SWATHS/PRS_L1_HCO/Data Fields/VNIR_Cube']
# vnir and swir are "class 'h5py._hl.Dataset'" meaning they are instances
# instance somehow inherits class characteristics in a class variable defined by the h5py module.
# In short, the SWIR array will have the same characteristics as the SWIR_Cube contained in immagine
# print(type(swir))
# print(vnir.shape)

# -------------------------------------- Find out the EPSG projection --------------------------------------
if tipo_prodotto=="PRS_L1_STD" or tipo_prodotto== "PRS_L2B_STD" or tipo_prodotto== "PRS_L2C_STD":
    # call the function to get the projection
    proiezione_epsg=proiezione(latitudine_pan,longitudine_pan)
    codice_epsg='non assegnato dal metadato'
    #proiezione_epsg=str(proiezione_epsg)
else:
    # the projection is reported in the metadata
    proiezione_epsg=int(codice_epsg)

proiez=('epsg:'+str(proiezione_epsg))

# ---------------------- Manually georeferencing if the variable 're-georeferencing' is 1 ----------------------
print ('If you have already created Geotiffs with this program')
print ('and you have verified a georeferencing error and')
print ('you have the new coordinates of the corners of the image then ')
ri_georeferenzia= int(input('insert 1 and then enter - otherwise press 2 or another number'))
if tipo_prodotto=="PRS_L2D_STD":
    sequenza = ['North_West NW ','South_West SW ','North_Est NE ','South_Est SE ']
    if ri_georeferenzia==1:
        print ('you need the geographic coordinate of the corners of the pan whole frame')
        print ("do\'nt look to the corners of the scene. Look to the frame")
else:
    sequenza = ['Top_Right NE ','Top_Left NW  ','Bottom_Right SE','Bottom_Left SW ']
    if ri_georeferenzia==1:
        print ('you need the geographic coordinate of the corners of the pan scene')
deltageoref=np.empty([6,4])
# deltageoref matrix is used to preserve the original values of the panchromatic vertices coordinates before modifying them.
# "ri_georeferenzia" should be left at 1 so that re_georeferencing the panchromatic will also re-georeference the hyperspectral. 
# For this reason the variables "ri_georeferenzia" and sequence are defined in the MAIN and not in the vertices function.

# -------------------------------------- window vertices table -------------------------------------- 
# call the vertex function to get the vertex table of the pan window and the size of the GSD
iper = 0
# assigns the variable iper = 0 because it is making the panchromatic
tabella_pan, x_pan_pixel, y_pan_pixel, deltageoref = vertici(latitudine_pan, longitudine_pan, proiezione_epsg, iper,tipo_prodotto,sequenza,ri_georeferenzia,deltageoref)
iper = 1
tabella_hyp, x_vnir_pixel, y_vnir_pixel, deltageoref = vertici(latitudine, longitudine, proiezione_epsg, iper,tipo_prodotto,sequenza,ri_georeferenzia,deltageoref)
# Please note in x_pan_Pixel it has the across dimension in the NO-NE direction, i.e. the E dimension
# Better would have been to call it e_pan_pixel and n_pan_pixel but that is how it is now

# ------------------------------------ writing the information file ------------------------------------
nome_file_tree= posizione_output + "informazioni"+data+".txt"
file_da_scrivere=open(nome_file_tree, 'w')
file_da_scrivere.write(nome_file_sorgente)
file_da_scrivere.write("\n Dati della finestra del PANcromatico \n")

file_da_scrivere.write("    "+ sequenza[0] +"    "+ sequenza[1] +"    "+ sequenza[2] +"    "+ sequenza[3]+"\n")
#if tipo_prodotto=="PRS_L1_STD" or tipo_prodotto== "PRS_L2B_STD" or tipo_prodotto== "PRS_L2C_STD":
#   file_da_scrivere.write("\n   Top_Right           Top_Left        Bottom_Right    Bottom_Left\n")
#else:
#    file_da_scrivere.write("\n   North_West         South_West       North_East      South_East\n")

for i in range (0,6):
    for j in range (0,4):
        if i <= 1:
            stringa = "      \t" + str(int(tabella_pan[i, j]))+"\t"
            file_da_scrivere.write(stringa)
        elif i==2 or i==3:
            file_da_scrivere.write("   \t")
            file_da_scrivere.write(format(tabella_pan[i, j],'.5f'))
        else:
            file_da_scrivere.write("   \t")
            file_da_scrivere.write(format(tabella_pan[i, j], '.2f'))
            #(format(tabella_pan[i, j], '.3f'))sta per 3 floatint("\n")
    file_da_scrivere.write('\n')

stringa=("\n"+"E GSD m ="+str(int(x_pan_pixel*10000000)/10000000)+ "  N GSD m ="+ str(int(y_pan_pixel*10000000)/10000000)+"\n")
file_da_scrivere.write(stringa)
file_da_scrivere.write("\n Dati dell'immagine cubica VNIR e SWIR \n")
file_da_scrivere.write("    "+ sequenza[0] +"    "+ sequenza[1] +"    "+ sequenza[2] +"    "+ sequenza[3]+"\n")
#if tipo_prodotto=="PRS_L1_STD" or tipo_prodotto== "PRS_L2B_STD" or tipo_prodotto== "PRS_L2C_STD":
#    file_da_scrivere.write("\n   Top_Right           Top_Left       Bottom_Right     Bottom_Left\n")
#else:
#    file_da_scrivere.write("\n   North_West         South_West       North_East      South_East\n")

for i in range (0,6):
    for j in range (0,4):
    #print (format(finestra[i, j], "10.2f"), "\t", end=" ")
        if i<=1:
            stringa = "    \t" + str(int(tabella_hyp[i, j]))+"\t"
            file_da_scrivere.write(stringa)
        elif i==2 or i==3:
            file_da_scrivere.write("   \t")
            file_da_scrivere.write(format(tabella_hyp[i, j],'.5f'))
        else:
            file_da_scrivere.write("   \t")
            file_da_scrivere.write(format(tabella_hyp[i, j],'.2f'))
    #print("\n")
    file_da_scrivere.write('\n')
stringa=("\n"+"E GSD m ="+str(int(x_vnir_pixel*1000000)/1000000)+ "  N GSD m ="+ str(int(y_vnir_pixel*1000000)/1000000)+"\n")
file_da_scrivere.write(stringa)
file_da_scrivere.close()

# ------------------------------------ Plot on screen pan and a couple of bands ------------------------------------

# the pyplot module of the matplot library (matplotlib.pyplot) has been renamed as plt
# the subplot command opens three parallel windows (1 line 3 windows)
# 14 inches wide and 6 inches high
# fig, axes = plt.subplots(1, 3, figsize=(17,7), sharex=True, sharey=True)
# sharex imposes the same axes. Not useful because one figure is 6000 and two are 1000
fig, axes = plt.subplots(1, 3, figsize=(17,7))
# Arguments: 1 is because it is figure 1, 3 because it is 3 subfigures.
# with "annotate" we output a text
# In xytext we tell it what percentage of the window from the bottom left. 
plt.annotate("Close this window!",xy=(.1,.1), xycoords='figure fraction', xytext=(.6,.9), fontsize=30, color='red')

# Start working on window axes[0].
plt.sca(axes[0])

# show panchromatic image, cmap=RdYlGn = red bottom green top,
# possible values for colourmap:
# https://matplotlib.org/2.0.1/examples/color/colormaps_reference.html
# massimo=np.max(pancromatico)
# mostra l'immagine pancromatico, cmap=RdYlGn = rosso in basso verde in alto,
minimo=np.min(pancromatico)
media=np.mean(pancromatico)
massimo=media*2

plt.imshow(pancromatico,cmap='Greys', vmin=minimo,vmax=massimo)
plt.colorbar(shrink=0.5)
plt.title("Panchromatic")
plt.xlabel('Column #')
plt.ylabel('Row #')


banda=47
imm_banda= estrai_banda(vnir,banda)
massimo=np.max(imm_banda)
minimo=np.min(imm_banda)
plt.sca(axes[1])
plt.imshow(imm_banda,cmap='Greens', vmin=minimo,vmax=massimo)
plt.colorbar(shrink=0.5)
plt.title("Green"+ str(list_cw_vnir[banda]) +" nm band "+ str(banda))
plt.xlabel('Column-along Track #')
plt.ylabel('Row #')

banda=26
imm_banda=estrai_banda(vnir,banda)
massimo=np.max(imm_banda)
minimo=np.min(imm_banda)
plt.sca(axes[2])
plt.imshow(imm_banda,cmap='Reds', vmin=minimo,vmax=massimo)
plt.colorbar(shrink=0.5)
plt.title("Red "+str( list_cw_vnir[banda] )+" nm band" + str(banda))
plt.xlabel('Column #')
plt.ylabel('Row #')
plt.show()
plt.close()

# ------------------------------------ create a small frame for the panchromatic ------------------------------------

# we create a layer consisting of a black area, as large as the panchromatic, with a thin slightly lighter frame around it
# element 0,1 in the table is 7283 and is the rows - height and we put it in J_max 
# in the table at index 1,2 there is 7355 the width of the image
j_max=int(tabella_pan[0, 1])
i_max=int(tabella_pan[1, 2])
# print ('imax  jmax,',i_max,j_max)
cornice=np.zeros ((j_max+1,i_max+1),dtype=np.uint16)
# in the example frame is an array of 0 of size 7284*7355
# the row is j_max=7283 the column is i_max 7355 but Python counts from 0 so you have to size 1 extra
# the next row is omitted
# val_corner=int(np.mean(panchromatic))
j=0
i=0
# in a pixel of the frame we assign a large value to display the colours of the frame
cornice[10,10]= 10000
for j in range(0,j_max+1):
    cornice[j,0]=4000
    # the first line will correspond to edge NO SO 
    cornice[j, i_max]=4350
    # the first line will correspond to edge NE SE 
for i in range (0,i_max+1):
    cornice[0,i]=4750
    # line will correspond to edge NO NE 
    cornice[j_max,i]=5000


# ------------------------ determination of the affinity matrix and its parameters ------------------------

''' L1 = data are already transformed into radiance and coregistered but not georeferenced. They are Top of atmosphere TOA.
L2b = obtained from L1 but atmospheric correction applied. They are radiance data. 
L2c = "at surface radiance" Angstrom correction, vapour correction and cloud thickness correction are applied. This is reflectance data.
L2d = geocoded geocoded i.e. DEM-corrected. '''

# In the case of L1/L2b/L2c data, the image must be rotated, mirrored and deformed.
# In the case of L2D data, the data are already reprojected and only need to generate the geotiff.
# we know the new NW coordinates but need the azimuth of the NW SW side (180° for L2D);
# in the table there are lists obtained from the data which are in the list table_hyp or table_pan;
# that is, two 4*6 matrices which for each row have the two indices, lat and long, and utm coordinates of the vertices; 
# in the lists of cardinal points below we enter pairs of UTM coordinates.
trasformazione_hyp=affinity(tipo_prodotto,tabella_hyp,x_vnir_pixel,y_vnir_pixel)

# ------------------------ Profile definition for affine transformation with Rasterio ------------------------

# Basically, in profile1 we define the parameters that we want to send to rasterio between the curly brackets.
profile1 = {
    'driver': 'GTiff',
    'dtype': 'uint16',
     #'nodata': 0,
    'width': vnir.shape[2],
    'height': vnir.shape[0],
    # width=Vnin.shape[2] gives the number of pixels of the along track dimension
    # the L1 level is rotated and reflected but essentially the width is across
    'tiled': True,
    # constructs the geotiff as a tiling and not by rows
    'count': 1,
    # count =1 is only one band
    'crs': {'init': proiez},
    'transform':trasformazione_hyp
}
# print ("profile1", profile1)

# ------------------------ writes the VNIR L1 bands projecting them with Rasterio ------------------------
# banda is 0 to 65 for vnir, 0 to 172 for swir and 0 for panchromatic
for banda in range (66):
    # crea il nome del file:
    tipo = "vnir"
    # "tipo" can assume: vnir or swir or pan
    lamda = (int(list_cw_vnir[banda] * 10))/10
    lung_onda=str(lamda)+"nm"
    nome_file_tree= posizione_output +"0"+lung_onda+codice+tipo+data+".tif"
    if banda==0:
        nome_file_tree = posizione_output + lung_onda +codice + tipo + data+".tif"
    print ("Now generating the GeoTIFF  ",nome_file_tree)
    print ("flag=",list_cw_vnir_flags[banda], "banda n", banda)
    # opens to write, using the affinity matrix parameters etc. specified in profile
    if list_cw_vnir_flags[banda]>0:
        with rasterio.open(nome_file_tree, "w",**profile1) as destinazione:
            destinazione.write(vnir[:, banda, :], indexes=1)
        destinazione.close()

# ------------------------ writes the SWIR L1 bands projecting them with rasterio ------------------------
for banda in range (173):
    tipo = "swir"
    lamda = (int(list_cw_swir[banda] * 10)) / 10
    lung_onda = str(lamda) + "nm"
    nome_file_tree= posizione_output +lung_onda+codice+tipo+data+".tif"
    if banda>=164:
        # adds a zero to the wavelength
        nome_file_tree = posizione_output + "0" + lung_onda+codice + tipo + data + ".tif"
    # opens to write, using the affinity matrix parameters etc. specified in profile
    print("flag=", list_cw_swir_flags[banda], "banda n", banda)
    # opens to write, using parameters specified in profile
    if list_cw_swir_flags[banda] > 0:
        print("Now generating the GeoTIFF  ", nome_file_tree)
        with rasterio.open(nome_file_tree, "w", **profile1) as destinazione:
            destinazione.write(swir[:, banda, :], indexes=1)
        destinazione.close()

# ------------------------ writes the PAN band projecting it with rasterio ------------------------
print ("Now writing the PAN GeoTIFF ")
trasformazione=affinity(tipo_prodotto,tabella_pan,x_pan_pixel,y_pan_pixel)
banda=0
tipo="pan"
profile2 = {
    'driver': 'GTiff',
    'dtype': 'uint16',
    #'nodata': 0,
    # 0-value data are NOT represented in the geotiff
    #'width': tabella_pan[0,1],
    #'height': tabella_pan[1,2],
    'width': pancromatico.shape[1],
    'height': pancromatico.shape[0],
    # note that the code projects the panchromatic image, pan[x,y] dimensions, into a [width,hight] window.
    # if we have dimensioned pan[x+1000,y+1000] larger, the program makes it all fit into the [width,hight]
    'count': 1,
    'tiled': True,
    # constructs the geotiff as tiled ('tiled') and not striped ('striped')
    'crs': {'init': proiez},
    #'crs': {'init': 'epsg:32633'},
    'transform': trasformazione
}
# print ("Transformation pan", trasformazione)
# print ("shapeW",pancromatico.shape[0],"shapeH",pancromatico.shape[1])
# print (x_pan_pixel,"x_pan_pixel",y_pan_pixel,'y_pan_pixel')
# print('width', tabella_pan[0,1],'height', tabella_pan[1,2])
lamda=0
nome_file_tree= posizione_output + "Pancromatico"+codice+data+".tif"
# opens to write, using parameters specified in profile
with rasterio.open(nome_file_tree, "w",**profile2) as destinazione:

    destinazione.write(pancromatico[:,:], indexes=1)
# print (profile2)

# ------------------------ Close file and exit ------------------------
destinazione.close()
sys.exit()