patrick::with_parameters_test_that(
  "Nonbatchwise WAveICA:",
  {
    input_data_path <- file.path("test-data", test_input)
    input_data <- readRDS(input_data_path)

    injection_order <- as.integer(input_data$injectionOrder)
    input_data <- dplyr::select(input_data, -any_of(c("sampleName", "injectionOrder", "sampleType", "batch", "class")))

    actual <- waveica_nonbatchwise(
      data = input_data,
      wf = "haar",
      injection_order = injection_order,
      K = 10,
      alpha = 0,
      cutoff = 0.1
    )

    expected_path <- file.path("test-data/nonbatchwise-correction", expected)
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
