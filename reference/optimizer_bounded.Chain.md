# Whether a Chain Takes Box Bounds

`TRUE` only when EVERY stage takes them: bounds are passed to all of
them, so one stage that would ignore them makes the chain unable to
promise the box.

## Arguments

- optimizer:

  A `Chain` object.

## Value

A single logical.
