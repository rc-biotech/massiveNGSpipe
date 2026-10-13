test_that('a partial reshift preserves unaffected offsets and keeps fresh replacements', {
  path<-tempfile(fileext='.rds');on.exit(unlink(path))
  old<-list('/data/A.ofst'=data.frame(fraction=28,offsets_start=-13),
            '/data/B.ofst'=data.frame(fraction=29,offsets_start=-12),
            '/data/excluded.ofst'=data.frame(fraction=30,offsets_start=-12))
  fresh<-list('/data/B.ofst'=data.frame(fraction=29,offsets_start=-14))
  saveRDS(fresh,path)
  pshift_preserve_sibling_offsets(path,old,c('/data/A.ofst'))
  result<-readRDS(path)
  expect_setequal(names(result),c('/data/A.ofst','/data/B.ofst'))
  expect_identical(result[['/data/A.ofst']],old[['/data/A.ofst']])
  expect_identical(result[['/data/B.ofst']],fresh[['/data/B.ofst']])
  expect_false('/data/excluded.ofst' %in% names(result))
  # A repeated repair must leave already-written fresh values intact.
  pshift_preserve_sibling_offsets(path,old,c('/data/A.ofst','/data/B.ofst'))
  expect_identical(readRDS(path),result)
})
