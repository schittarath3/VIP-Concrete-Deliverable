<h1> VIP Concrete Image Processing </h1>

Overview: Repoistory for code processing binary images into .stl files of concrete

<h2> Repo Contents: </h2>

1. Image Processing Folder: Code to process CT scans of concrete into .stl files.
    1. image_splitter.py - Python code used to split each original image into the four quadrants
        1. Requirements:
            1. Python 3.7.6
            2. opencv2
            3. os
            4. numpy

    2. split_to_binary - MATLAB code to convert split images into binary images of aggregates
        1. Requirements
            1. MATLAB (version used was 2019b)
            2. Image Batch Processor App (part of Image Processing and Computer Vision package)
    3. Other programs:
        + Fiji- ImageJ distribution with many plugins for scientific analysis
        + (Optional) Meshmixer

2. Main Folder: Current working directory
    Current pipeline:
        1. reducemesh.m
        2. cubelets.m
        3. insertAgg.m
        4. growAgg.m
        5. aggCubeTest.m

Note: Detailed description of the functions of each code is described in comments
