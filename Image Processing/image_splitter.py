#load dependencies
import cv2
import os
import numpy as np


def get_contours(species_image):
        image = cv2.imread(species_image)
        gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
        inverted = (255-gray)

        #create threshold; original min value = 215
        thresh = cv2.threshold(inverted, 195,255, cv2.THRESH_TOZERO_INV)[1]

        # Find contours
        cnts = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)[0]
        return cnts

def image_split(image, mainPath, image_name, cnts):
    image = cv2.imread(image)
    original = image.copy()
    
    #create output directories
    outpath = os.path.join(mainPath, "Final Images")
    outpath_bl = os.path.join(outpath, "bottom_left")
    outpath_br = os.path.join(outpath, "bottom_right")
    outpath_tl = os.path.join(outpath, "top_left")
    outpath_tr = os.path.join(outpath, "top_right")
    if not os.path.exists(outpath):
        os.makedirs(outpath)
        os.makedirs(outpath_bl)
        os.makedirs(outpath_br)
        os.makedirs(outpath_tl)
        os.makedirs(outpath_tr)

    #Iterate thorugh contours and returns split images
    thold_area = 5000 #area of image quadrants
    image_number = 1
    for c in cnts:
        area = cv2.contourArea(c)         
        if area > thold_area: #filters out small objects
            rect = cv2.minAreaRect(c)
            box = np.int0(cv2.boxPoints(rect))
            width = int(rect[1][0])
            height = int(rect[1][1])
            src_pts = box.astype("float32")
            dst_pts = np.array([[0, height-1],
                            [0, 0],
                            [width-1, 0],
                            [width-1, height-1]], dtype="float32")
            matrix = cv2.getPerspectiveTransform(src_pts, dst_pts)
            warped = cv2.warpPerspective(original, matrix, (width, height))
            if image_number  == 1:
                cv2.imwrite((os.path.join(outpath_bl, "bl_{}".format(image_name))), warped)
            elif image_number == 2:
                cv2.imwrite((os.path.join(outpath_br, "br_{}".format(image_name))), warped)
            elif image_number == 3:
                cv2.imwrite((os.path.join(outpath_tl, "tl_{}".format(image_name))), warped)
            elif image_number == 4:
                cv2.imwrite((os.path.join(outpath_tr, "tr_{}".format(image_name))), warped)
                image_number = 1
            image_number += 1
        
#main function
def master(inPath, mainPath, species_image):
    cnts = get_contours(species_image)
    for image_name in os.listdir(inPath):  #reads all images in input folder, inPath
        input_path = os.path.join(inPath, image_name)
        image_split(input_path, mainPath, image_name, cnts)

master(r"path\to\folder", r"path\to\output\folder", r"path\to\ideal\image")
