# Uncollect files

Given a file with many different slots, create multipe different files with a dedicated subset of slots.
This is particularly useful when you want to undo the collect module and compute downstream task on different slots (layers, obsm, etc) separately.

## Parameters

```yaml
DATASETS:
  test:
    input:
      uncollect: test/input/
    uncollect:
      new_file_ids:
        - file_1
        - file_2
        - file_3
      sep: --
```

* `sep`: separator used for matching `new_file_id` in the slots. Note, the separator should be unique to the slots you want to match, to avoid unintended side effects. Default: '--'
* `new_file_ids`: list of new file ids to create. If the file name (`new_file_id`) uniquely matches a slot name, that slot will be included in the new file. For slots that don't contain the unique separator, they will be included in **all new files**. Only in cases where the separator is in the slot name, but the slot does not match the `new_file_id`, the slot will not be included.
