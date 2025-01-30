// Prompt user to open an image
showMessage("Select an image for analysis");
open();

// Prompt user to measure the scale bar manually using the Magic Wand
showMessage("Zoom into scale bar and use the Magic Wand tool to measure.");

// Wait for the user to perform the measurement
waitForUser("Once you've measured the scale bar, click OK.");

// Retrieve the "Major" length from the last result
measuredLength = getResult("Major", nResults-1);

// Prompt user for the known length of the scale bar
Dialog.create("Input Scale Bar Information");
Dialog.addChoice("Scale bar type:", newArray("1 mm", "500 µm"), "1 mm");
Dialog.show();

// Get user input
scaleBarType = Dialog.getChoice();

// Determine known length based on scale bar type
if (scaleBarType == "1 mm") {
    knownLength = 1; // in mm
} else if (scaleBarType == "500 µm") {
    knownLength = 1; // Scale equivalent to 1 mm
    measuredLength = measuredLength * 2; // Adjust measurement to match 1 mm scale
}

// Automatically set the scale globally
run("Set Scale...", "distance=" + measuredLength + " known=" + knownLength + " unit=mm");

// Check set scale box
run("Set Scale...");

// Deselect the current selection
run("Select None");

// Zoom out the image
run("Scale to Fit");

// Duplicate the image
run("Duplicate...", "title=Duplicate");

// Black and white original image conversion
waitForUser("Select the original image, then click OK to proceed.");
run("RGB Stack");
run("Rotate 90 Degrees Left");


// Split channels
waitForUser("Select the duplicated image, then click OK to split channels.");
run("Split Channels");

// Rename and assign channels
waitForUser("Select the DAPI (blue) channel window, then click OK.");
rename("DAPI");
waitForUser("Select the GFP (green) channel window, then click OK.");
rename("GFP");
waitForUser("Select the mCherry (red) channel window, then click OK.");
rename("mCherry");

// Merge channels
run("Merge Channels...", "c1=mCherry c2=GFP c3=DAPI create");
run("Rotate 90 Degrees Left");


// Open channels tool for user interaction
run("Channels Tool...");
waitForUser("Use the Channels Tool to adjust select red and green channels, then click OK to finish.");