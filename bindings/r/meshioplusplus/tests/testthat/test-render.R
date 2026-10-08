# Software rendering (v16.33.0, field rendering v16.34.0).

cube_surface <- function() {
  pts <- matrix(c(
    0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
    0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1
  ), nrow = 3)
  quads <- matrix(c(
    1, 4, 3, 2, 5, 6, 7, 8, 1, 2, 6, 5,
    4, 8, 7, 3, 1, 5, 8, 4, 2, 3, 7, 6
  ), nrow = 4)
  m <- mio_mesh()
  mio_set_points(m, pts)
  mio_add_cell_block(m, "quad", quads)
  mio_add_point_data(m, "x", pts[1, ])
  m
}

test_that("mio_render returns a frame with pixels, ids and notes", {
  m <- cube_surface()
  on.exit(mio_release(m))
  f <- mio_render(m, width = 48L, height = 32L, background = mio_rgba(255, 255, 255))
  expect_s3_class(f, "mio_frame")
  expect_equal(f$width, 48L)
  expect_equal(f$height, 32L)
  expect_equal(dim(f$rgba), c(4L, 48L, 32L))
  expect_equal(dim(f$cell_ids), c(48L, 32L))
  expect_equal(as.integer(f$rgba[, 1, 1]), c(255L, 255L, 255L, 255L))
  expect_gt(sum(f$cell_ids >= 0), 100)
  expect_null(f$range)
  expect_equal(dim(mio_frame_array(f)), c(32L, 48L, 4L))
  expect_equal(mio_frame_array(f)[1, 1, 1], 1)
})

test_that("a mapped field reports its range, ticks and symmetric limits", {
  m <- cube_surface()
  on.exit(mio_release(m))
  g <- mio_render(m,
    width = 40L, height = 40L, color_by = "x", colorbar = TRUE, symmetric = TRUE,
    clip = c(5, 95), view = "+z"
  )
  expect_false(is.null(g$range))
  expect_equal(g$range[1], -g$range[2])
  expect_true(any(startsWith(g$notes, "x: ")))
  expect_true(any(startsWith(g$notes, "ticks: ")))
})

test_that("text, PNG and file forms render", {
  m <- cube_surface()
  on.exit(mio_release(m))
  txt <- mio_render_text(m, cols = 20L, rows = 6L, color_depth = "mono")
  expect_equal(lengths(regmatches(txt, gregexpr("\n", txt))), 6L)
  expect_false(grepl("\033", txt, fixed = TRUE))
  expect_true(startsWith(mio_render_text(m, cols = 20L, rows = 6L, format = "html"), "<!DOCTYPE html>"))

  png <- mio_render_png(m, width = 30L, height = 20L)
  expect_type(png, "raw")
  expect_equal(png[1:8], as.raw(c(0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a)))
  expect_error(mio_render_png(m, compress = 10L))

  path <- tempfile(fileext = ".png")
  mio_write_snapshot(path, m, width = 32L, height = 24L)
  expect_equal(readBin(path, "raw", 4), as.raw(c(0x89, 0x50, 0x4e, 0x47)))
  expect_error(mio_write_snapshot(tempfile(fileext = ".vtu"), m))
})

test_that("bad options are refused by name", {
  m <- cube_surface()
  on.exit(mio_release(m))
  expect_error(mio_render(m, color_by = "nope"), "nope")
  expect_error(mio_render(m, bogus = 1), "unknown option 'bogus'")
  expect_error(mio_render(m, shading = "plastic"), "shading")
  expect_error(mio_render(m, width = 0L), "width and height")
})

test_that("the terminal colour depth follows the environment variables", {
  expect_equal(mio_detect_color_depth(no_color = "1"), "mono")
  expect_equal(mio_detect_color_depth(color_term = "truecolor"), "truecolor")
  expect_equal(mio_detect_color_depth(term = "xterm-256color"), "256")
  expect_equal(mio_detect_color_depth(term = "xterm"), "16")
  expect_equal(mio_rgba(200, 197, 189), 0xC8C5BDFF)
})
