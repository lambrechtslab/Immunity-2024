import qupath.lib.gui.dialogs.Dialogs
def pathModel2 = Dialogs.showMessageDialog("Please select your stardist segmentation file","Please select your StarDist cell segmentation file")
def pathModel = Dialogs.promptForFile(null)
def pathModel3 = pathModel.toString()
//def pathModel = Dialogs.showInputDialog("My title", "Please paste the StarDist file path without quotes", 'text')
import qupath.ext.stardist.StarDist2D

// if you are always using the same StarDist model you can set the path below instead of selecting via dialog box
//def pathModel = "C:\\Users\\gcarlson\\OneDrive - Akoya Biosciences\\Desktop\\CODEX Analysis Demo\\CODEX_cell_seg.pb"

def stardist = StarDist2D.builder(pathModel3)
        .threshold(0.5)              // Probability (detection) threshold
        .channels(0)            // Select detection channel
        .normalizePercentiles(1, 99) // Percentile normalization
        .pixelSize(0.325)              // Resolution for detection
        .cellExpansion(5.0)          // Approximate cells based upon nucleus expansion
        .cellConstrainScale(1.5)     // Constrain cell expansion using nucleus size
        .measureShape()              // Add shape measurements
        .measureIntensity()          // Add cell measurements (in all compartments)
        .includeProbability(true)    // Add probability as a measurement (enables later filtering)
        .build()

// Run detection for the selected objects
def imageData = getCurrentImageData()
def pathObjects = getSelectedObjects()
if (pathObjects.isEmpty()) {
    Dialogs.showErrorMessage("StarDist", "Please select a parent object!")
    return
}
stardist.detectObjects(imageData, pathObjects)
println 'Done!'