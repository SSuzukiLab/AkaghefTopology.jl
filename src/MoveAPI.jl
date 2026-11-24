"""Generic names for moves and for the search database. Presentations add methods; this module has no implementation."""
module MoveAPI

"""
    movable(x)
    movable(x, kind; choices...)

Witnesses for the moves that apply to `x`. `movable(x)` returns every witness
the presentation currently offers, and that collection is finite.
`movable(x, kind; choices...)` keeps one kind. An omitted choice means every
candidate of that kind.
"""
function movable end

"""
    move(x, witness)

Return the value of applying `witness` to `x`. `x` is left unchanged.
"""
function move end

"""
    move!(x, ...)

Destructive `move`. Only a mutable presentation (`VirtualLink`) defines it.
"""
function move! end

"""
    search_state(x)

Key under which `x` is stored in a search database. Equal presentations share a key.
"""
function search_state end

"Store a search-database state."
function put_search_state! end

"Load a stored search-database state."
function load_search_state end

"Record the start of a search run."
function begin_search_run! end

"Record the end of a search run."
function finish_search_run! end

"Apply the selected moves of `x` and store the resulting transitions."
function expand_search_state! end

"Read the stored neighbors of a search state."
function search_neighbors end

"Read stored paths in a search database."
function query_search_paths end

"Replay every stored transition and check the search database."
function verify_search_store end

export movable, move, move!,
       search_state, put_search_state!, load_search_state,
       begin_search_run!, finish_search_run!, expand_search_state!,
       search_neighbors, query_search_paths, verify_search_store

end
