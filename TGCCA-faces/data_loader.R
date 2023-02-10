data_loader <- function(subject, pose = 51, illumination = 2, expression = 1,
                        size = 64) {
  subject <- sprintf("%03d", subject)
  pose <- sprintf("%03d", pose)
  illumination <- sprintf("%02d", illumination)
  expression <- sprintf("%02d", expression)

  data_path <- paste0(
    "./data/", subject, "/", subject, "_01_", expression,
    "_", pose, "_", illumination, "_crop_128.png"
  )

  if (file.exists(data_path)) {
    im <- imager::grayscale(imager::load.image(data_path))
    return(imager::resize(im, size, size, interpolation_type = 3))
  } else {
    stop(
      "File ", data_path, " does not exist. Did you download the dataset and ",
      "put the images in the data folder?"
    )
  }
}
