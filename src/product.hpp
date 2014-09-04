// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#ifndef product_hpp
#define product_hpp

namespace flexiblesusy {

#define PRODUCT(...) (get_product(__VA_ARGS__)(__VA_ARGS__))

#define get_product(...) \
    get_product_macro(__VA_ARGS__, product_user_t, product_size_t,)

#define get_product_macro(_1, _2, _3, _4, _5, name, ...) name

#define product_size_t(idx, ini, fin, expr)		\
    product<size_t, (ini), (fin)>([&](size_t (idx)) { return (expr); })

#define product_user_t(type, idx, ini, fin, expr)	\
    product<type, (ini), (fin)>([&](type (idx)) { return (expr); })

template<class Idx, Idx ini, Idx fin, class Function>
auto product(Function f) -> decltype(f(ini))
{
    decltype(f(ini)) s = 1;
    for (Idx i = ini; i <= fin; i++) s *= f(i);
    return s;
}

#define UPRODUCT(...) (get_uproduct(__VA_ARGS__)(__VA_ARGS__))

#define get_uproduct(...) \
    get_uproduct_macro(__VA_ARGS__, uproduct_user_t, uproduct_size_t,)

#define get_uproduct_macro(_1, _2, _3, _4, _5, name, ...) name

#define uproduct_size_t(idx, ini, fin, expr)		\
    uproduct<size_t, (ini), (fin)>([&](size_t (idx)) { return (expr); })

#define uproduct_user_t(type, idx, ini, fin, expr)	\
    uproduct<type, (ini), (fin)>([&](type (idx)) { return (expr); })

template<bool valid, class Function, class Idx, Idx ini, Idx fin>
struct unroll_product;

template<class Function, class Idx, Idx ini, Idx fin>
struct unroll_product<false, Function, Idx, ini, fin> {
    static auto eval(Function f) -> decltype(f(ini)) {
	return decltype(f(ini))(1);
    }
};

template<class Function, class Idx, Idx idx>
struct unroll_product<true, Function, Idx, idx, idx> {
    static auto eval(Function f) -> decltype(f(idx)) {
	return f(idx);
    }
};

template<class Function, class Idx, Idx ini, Idx fin>
struct unroll_product<true, Function, Idx, ini, fin> {
    static const Idx mid = (ini+fin)/2;
    static auto eval(Function f) -> decltype(f(mid)) {
	return
	    unroll_product<(mid > ini), Function, Idx, ini, mid-1>::eval(f) *
	    f(mid) *
	    unroll_product<(fin > mid), Function, Idx, mid+1, fin>::eval(f);
    }
};

template<class Idx, Idx ini, Idx fin, class Function>
auto uproduct(Function f) -> decltype(f(ini))
{
    return unroll_product<(fin >= ini), Function, Idx, ini, fin>::eval(f);
}

} // namespace flexiblesusy

#endif // product_hpp
