#!/bin/bash
mylist="PDBDEV_00000001_resize, PDBDEV_00000002, PDBDEV_00000003, PDBDEV_00000004, PDBDEV_00000005, PDBDEV_00000006, PDBDEV_00000007, PDBDEV_00000008, PDBDEV_00000009, PDBDEV_00000010, PDBDEV_00000011, PDBDEV_00000012, PDBDEV_00000014, PDBDEV_00000015, PDBDEV_00000016, PDBDEV_00000017, PDBDEV_00000018, PDBDEV_00000020, PDBDEV_00000021, PDBDEV_00000022, PDBDEV_00000023, PDBDEV_00000024, PDBDEV_00000025, PDBDEV_00000026, PDBDEV_00000027, PDBDEV_00000028, PDBDEV_00000029, PDBDEV_00000031, PDBDEV_00000032, PDBDEV_00000033"
Field_Separator=$IFS

IFS=", "
for filename in $mylist; do
    echo $filename
    file=$filename + ".png"
    echo $file
    #convert $filename  +profile "*" -fuzz 1% -trim +repage -resize 400x400 -background white -gravity center -extent 400x400 $filename_resize.png

done
