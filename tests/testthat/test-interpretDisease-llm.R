test_that("llm explanation serializes structured evidence and records audit metadata", {
    x <- interpretDisease(c("1", "2", "9", "10"))
    captured_request <- NULL

    adapter <- function(request) {
        captured_request <<- request
        list(
            text = "Supported by the supplied overlap evidence.",
            cited_evidence_ids = request$payload$evidence$evidence_id[[1]],
            unsupported_claims = character(),
            metadata = list(provider = "fake", finish_reason = "stop")
        )
    }

    x <- explainDisease(
        x,
        method = "llm",
        provider = "fake-provider",
        model = "fake-model",
        language = "en",
        style = "research",
        generation_params = list(temperature = 0),
        llm_adapter = adapter
    )

    expect_equal(x@explanation$method, "llm")
    expect_length(x@explanation$text, 1)
    expect_true(grepl("You may only use the evidence supplied below.", captured_request$prompt, fixed = TRUE))
    expect_equal(captured_request$payload$query$query_id[[1]], "q1")
    expect_true(nrow(captured_request$payload$evidence) > 0)

    metadata <- x@explanation$metadata$queries$q1
    expect_equal(x@explanation$metadata$provider, "fake-provider")
    expect_equal(x@explanation$metadata$model, "fake-model")
    expect_equal(x@explanation$metadata$language, "en")
    expect_equal(x@explanation$metadata$style, "research")
    expect_equal(metadata$generation_params$temperature, 0)
    expect_true(nzchar(metadata$payload_hash))
    expect_true(nzchar(metadata$prompt_hash))
    expect_true(nzchar(metadata$response_hash))
    expect_true(nzchar(metadata$evidence_hash))
    expect_equal(metadata$cited_evidence_ids[[1]], captured_request$payload$evidence$evidence_id[[1]])
})

test_that("llm explanation rejects unsupported evidence claims", {
    x <- interpretDisease(c("1", "2", "9", "10"))

    expect_error(
        explainDisease(
            x,
            method = "llm",
            llm_adapter = function(request) {
                list(
                    text = "This adds an unsupported fact.",
                    cited_evidence_ids = request$payload$evidence$evidence_id[[1]],
                    unsupported_claims = "external unsupported fact"
                )
            }
        ),
        "unsupported claims"
    )
})

test_that("llm explanation surfaces adapter failures clearly", {
    x <- interpretDisease(c("1", "2", "9", "10"))

    expect_error(
        explainDisease(
            x,
            method = "llm",
            llm_adapter = function(request) {
                stop("provider timeout")
            }
        ),
        "LLM adapter failed"
    )
})

test_that("interpretDisease can call the llm explanation convenience path", {
    x <- interpretDisease(
        c("1", "2", "9", "10"),
        explain = "llm",
        llm_adapter = function(request) {
            list(
                text = "Convenience path explanation.",
                cited_evidence_ids = request$payload$evidence$evidence_id[[1]],
                unsupported_claims = character(),
                metadata = list(path = "interpretDisease")
            )
        }
    )

    expect_equal(x@explanation$method, "llm")
    expect_match(x@explanation$text[[1]], "Convenience path explanation.", fixed = TRUE)
    expect_equal(x@explanation$metadata$queries$q1$response_metadata$path, "interpretDisease")
})
