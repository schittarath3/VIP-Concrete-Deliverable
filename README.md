<h1> VIP Concrete Image Processing </h1>

Overview: Code to process CT scans of concrete into .stl files.

<h2> Repo Contents: </h2>

1. Example Folder: Contains input images and outputs
	1. Original Images- Contain original images given from VIP Mentor
	2. Split Images- output from image_splitter.py
	3. Binary Images - output from split_to_binary.m
	4. bl_54_94.stl - output from Fiji (ImageJ) .stl creator; encoded in binary
	5. bl_54_94_MM.stl - output from Meshmixer after removing artifacts and reducing
	
2. image_splitter.py - Python code used to split each original image into the four quadrants
	1. Requirements:
		1. Python 3.7.6
		2. opencv2
		3. os
		4. numpy
		
3. split_to_binary - MATLAB code to convert split images into binary images of aggregates
	1. Requirements
		1. MATLAB (version used was 2019b)
		2. Image Batch Processor App (part of Image Processing and Computer Vision package)

Other programs:
+ Fiji- ImageJ distribution with many plugins for scientific analysis
+ (Optional) Meshmixer

Note: Detailed description of the functions of each code is described in comments
