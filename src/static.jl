# Minor functionality from Static.jl, reimplemented here to avoid invalidations
struct _StaticInt{N} end
_StaticInt(N::Int) = _StaticInt{N}()
const _IntOrStaticInt = Union{Integer, _StaticInt}
_dim(::Type{_StaticInt{N}}) where {N} = N::Int
Base.:(==)(a::Number, b::_StaticInt{N}) where {N} = a == N
Base.:(==)(a::_StaticInt{N}, b::Number) where {N} = b == N

struct _True end
struct _False end
_static(b::Bool) = b ? _True() : _False()