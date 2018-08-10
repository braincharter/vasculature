filename = getArgument();

input = filename + ".nii";
print("Reading image : " + input);
open(input);

print("Convert image in 8-bit");
run("8-bit");

print("Run Skel 3D");
run("Skeletonise 3D");

output = filename + "_skel.nii";
print("Save image : " + output);
run("NIfTI-1", "save=&output");

exit("Completed")
