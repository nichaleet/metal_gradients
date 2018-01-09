pro maketif,name,red,green,blue
red = readfits(red)
green = readfits(green)
blue = readfits(blue)

write_tiff,name,red=Red, GREEN=Green, BLUE=Blue,PLANARCONFIG=2

end
