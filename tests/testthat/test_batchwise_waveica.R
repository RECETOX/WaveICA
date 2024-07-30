patrick::with_parameters_test_that(
  "Batchwise WaveICA:",
  {
    input_data_path <- file.path("test-data", test_input)
    input_data <- readRDS(input_data_path)

    batch <- as.integer(input_data$batch)
    group <- as.integer(input_data$sampleType)
    input_data <- dplyr::select(input_data, -any_of(c("sampleName", "injectionOrder", "sampleType", "batch", "class")))

    actual <- waveica(
      data = input_data,
      wf = "haar",
      group = group,
      batch = batch,
      K = 20,
      batch_threshold = 0.05,
      group_threshold = 0.05,
      alpha = 0
    )

    expected_path <- file.path("test-data/batchwise-correction", expected)
    expected <- readRDS(expected_path)

    expect_equal(actual, expected)
  },
  patrick::cases(
    amide_v1 = list(
      test_input = "amide_data_v1.rds",
      expected = "amide_corrected_v1.rds"
    ),
    amide_v2 = list(
      test_input = "amide_data_v2.rds",
      expected = "amide_corrected_v2.rds"
    )
  )
)
