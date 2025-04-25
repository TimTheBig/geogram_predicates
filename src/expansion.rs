use smallvec::SmallVec;
use core::{cmp::Ordering, fmt};

#[derive(Clone)]
pub struct Expansion {
    /// capacity is dynamic but starts inline
    data: SmallVec<[f64; 2]>,
}

/// Build an `Expansion` from a comma-separated list of `f64` literals or expressions.
///
/// ```rust
/// # use your_crate::Expansion;
/// let e = expansion![1.0, 2.5, 3.75];
/// assert_eq!(e.length(), 3);
/// assert_eq!(&e[..], &[1.0, 2.5, 3.75]);
/// ```
#[macro_export]
macro_rules! expansion {
    // match zero or more comma-separated expressions, allow trailing comma
    ($($x:expr),* $(,)?) => {
        // create a fixed-length array `[f64; count]`, then call From<[f64;N]>
        $crate::Expansion::from([$($x),*])
    };
    ($x:expr) => {
        $crate::Expansion::from($x)
    };
}

impl Default for Expansion {
    fn default() -> Self {
        Self { data: SmallVec::new() }
    }
}

impl Expansion {
    /// Create a new `Expansion` with the given capacity.
    ///
    /// Internally uses a `SmallVec<[f64;2]>` so small expansions
    /// are stored inline; larger ones spill to the heap.
    ///
    /// # Parameters
    /// - `capacity`: maximum number of components (not including
    ///   the extra sentry slot).
    ///
    /// # Examples
    /// ```
    /// let e = Expansion::with_capacity(5);
    /// assert_eq!(e.capacity(), 5);
    /// assert_eq!(e.length(), 0);
    /// ```
    pub fn with_capacity(capacity: usize) -> Self {
        let mut data = SmallVec::with_capacity(capacity + 1);
        data.resize(capacity + 1, 0.0);
        Self { data }
    }

    pub fn len(&self) -> usize {
        self.data.len()
    }

    pub fn capacity(&self) -> usize {
        self.data.capacity()
    }

    pub fn data(&self) -> &[f64] {
        &self.data
    }

    pub const fn data_mut(&mut self) -> &mut SmallVec<[f64; 2]> {
        &mut self.data
    }

    pub fn get(&self, i: usize) -> f64 {
        debug_assert!(i < self.data.len());
        self.data[i]
    }

    pub fn get_mut(&mut self, i: usize) -> &mut f64 {
        debug_assert!(i < self.data.len());
        &mut self.data[i]
    }

    /// Assign a single double to this expansion.
    ///
    /// After this call, `self.length() == 1` and `self[0] == a`.
    ///
    /// # Examples
    /// ```
    /// let mut e = Expansion::with_capacity(2);
    /// e.assign(3.14);
    /// assert_eq!(e.length(), 1);
    /// assert_eq!(e[0], 3.14);
    /// ```
    pub fn assign(&mut self, a: f64) -> &mut Self {
        self.data.clear();
        self.data.push(a);
        self
    }

    pub fn assign_expansion(&mut self, other: &mut Expansion) -> &mut Self {
        self.data.append(other.data_mut());
        for i in 0..other.len() {
            self.data[i] = other.data[i];
        }
        self
    }

    pub fn assign_abs(&mut self, rhs: &mut Expansion) -> &mut Self {
        self.assign_expansion(rhs);
        if self.sign() == Sign::Negative {
            self.negate();
        }
        self
    }

    /// Negate every component of the expansion in place.
    ///
    /// # Examples
    /// ```
    /// let mut e = Expansion::from(2.0);
    /// let f = -e.clone();
    /// assert_eq!(f[0], -2.0);
    /// ```
    pub fn negate(&mut self) -> &mut Self {
        for v in self.data.iter_mut() {
            *v = -*v;
        }
        self
    }

    pub fn scale_fast(&mut self, s: f64) -> &mut Self {
        for v in self.data.iter_mut() {
            *v *= s;
        }
        self
    }

    /// Estimate the value of this expansion by summing all components.
    ///
    /// This gives a quick—and not fully accurate—“approximate” value.
    ///
    /// # Examples
    /// ```
    /// let mut e = Expansion::with_capacity(3);
    /// e.assign(1.0);
    /// e.data_mut()[1] = 0.0000001;
    /// assert!(e.estimate() > 1.0);
    /// ```
    pub fn estimate(&self) -> f64 {
        self.data.iter().sum()
    }

    pub fn sign(&self) -> Sign {
        if self.len() == 0 {
            Sign::Zero
        } else {
            geo_sgn(*self.data.last().unwrap())
        }
    }

    pub fn equals(&self, rhs: &Expansion) -> bool {
        self.compare(rhs) == Sign::Zero
    }

    /// Compare two expansions by their estimated values.
    ///
    /// Returns:
    /// - `Ordering::Less` if `self.estimate() < other.estimate()`,
    /// - `Ordering::Equal` if they are (approximately) equal,
    /// - `Ordering::Greater` otherwise.
    ///
    /// # Examples
    /// ```
    /// let a = Expansion::from(1.0);
    /// let b = Expansion::from(2.0);
    /// assert!(a < b);
    /// ```
    pub fn compare(&self, rhs: &Expansion) -> Sign {
        let est_self = self.estimate();
        let est_rhs = rhs.estimate();
        geo_sgn(est_self - est_rhs)
    }
}

// Helper for sign handling
#[derive(PartialEq, Debug)]
pub(crate) enum Sign {
    Negative,
    Zero,
    Positive,
}

impl From<Sign> for i8 {
    fn from(sign: Sign) -> i8 {
        match sign {
            Sign::Positive => 1,
            Sign::Zero => 0,
            Sign::Negative => -1,
        }
    }
}

fn geo_sgn(value: f64) -> Sign {
    if value > 0.0 {
        Sign::Positive
    } else if value < 0.0 {
        Sign::Negative
    } else {
        Sign::Zero
    }
}

impl fmt::Debug for Expansion {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Expansion[{}] = [", self.len())?;
        for x in &self.data {
            write!(f, "{} ", x)?;
        }
        write!(f, "]")
    }
}

impl PartialEq for Expansion {
    fn eq(&self, other: &Self) -> bool {
        self.data == other.data
    }
}

impl PartialOrd for Expansion {
    fn partial_cmp(&self, other: &Expansion) -> Option<Ordering> {
        // Always returns Some(Ordering) because compare() gives a total order.
        let est_self = self.estimate();
        let est_rhs = other.estimate();
        est_self.partial_cmp(&est_rhs)
    }
}

impl core::ops::Index<usize> for Expansion {
    type Output = f64;
    fn index(&self, idx: usize) -> &f64 {
        &self.data[idx]
    }
}

impl From<f64> for Expansion {
    /// Create a length-1 expansion holding exactly `a`.
    fn from(a: f64) -> Self {
        let mut e = Expansion::with_capacity(1);
        e.assign(a);
        e
    }
}

impl<const N: usize> From<[f64; N]> for Expansion {
    /// Create an expansion from an array of `N` doubles.
    ///
    /// The resulting expansion has length `N` and
    /// components equal to the array elements in order.
    fn from(arr: [f64; N]) -> Self {
        let mut e = Expansion::with_capacity(N);
        // push each element in turn
        for &v in &arr {
            // we know capacity ≥ N, so this cannot panic
            e.data_mut().push(v);
        }
        e
    }
}

impl From<&[f64]> for Expansion {
    /// Create an expansion from a slice of doubles.
    ///
    /// The resulting expansion has length `slice.len()` and
    /// components equal to the slice elements in order.
    fn from(slice: &[f64]) -> Self {
        Expansion { data: SmallVec::from_slice(slice) }
    }
}

/// Compute the 3×3 determinant of nine `Expansion`s, returning a new `Expansion`.
///  
/// This mirrors the C++ macro:
/// `expansion_det3x3(a11, …, a33)` →
/// `new_expansion_on_stack(det3x3_capacity(...))->assign_det3x3(...)`  
#[macro_export]
macro_rules! expansion_det3x3 {
    (
        $a11:expr, $a12:expr, $a13:expr,
        $a21:expr, $a22:expr, $a23:expr,
        $a31:expr, $a32:expr, $a33:expr
    ) => {{
        // Compute exactly the capacity needed
        let cap = $crate::Expansion::det3x3_capacity(
            &$a11, &$a12, &$a13,
            &$a21, &$a22, &$a23,
            &$a31, &$a32, &$a33
        );
        // Allocate an Expansion with that capacity
        let mut e = $crate::Expansion::with_capacity(cap);
        // Perform the determinant assignment
        e.assign_det3x3(
            &$a11, &$a12, &$a13,
            &$a21, &$a22, &$a23,
            &$a31, &$a32, &$a33
        );
        e
    }};
}

impl Expansion {
    /// Compute the capacity needed to form the 3×3 determinant
    /// of the nine expansions a11…a33.
    ///
    /// Mirrors C++:
    /// ```
    /// index_t c11 = det2x2_capacity(a22,a23,a32,a33);
    /// index_t c12 = det2x2_capacity(a21,a23,a31,a33);
    /// index_t c13 = det2x2_capacity(a21,a22,a31,a32);
    /// return 2*(a11.len()*c11 + a12.len()*c12 + a13.len()*c13);
    /// ```
    pub fn det3x3_capacity(
        a11: &Expansion, a12: &Expansion, a13: &Expansion,
        a21: &Expansion, a22: &Expansion, a23: &Expansion,
        a31: &Expansion, a32: &Expansion, a33: &Expansion
    ) -> usize {
        let c11 = Self::det2x2_capacity(a22, a23, a32, a33);
        let c12 = Self::det2x2_capacity(a21, a23, a31, a33);
        let c13 = Self::det2x2_capacity(a21, a22, a31, a32);
        2 * (
            a11.length() * c11 +
            a12.length() * c12 +
            a13.length() * c13
        )
    }

    /// Assign to `self` the 3×3 determinant of the nine expansions.
    ///
    /// Computes
    /// ```txt
    /// a11·det2x2(a22,a23,a32,a33)
    /// − a12·det2x2(a21,a23,a31,a33)
    /// + a13·det2x2(a21,a22,a31,a32)
    /// ```
    pub fn assign_det3x3(
        &mut self,
        a11: &Expansion, a12: &Expansion, a13: &Expansion,
        a21: &Expansion, a22: &Expansion, a23: &Expansion,
        a31: &Expansion, a32: &Expansion, a33: &Expansion
    ) -> &mut Self {
        // 1) build the three 2×2 minors
        let mut m11 = Expansion::with_capacity(
            Self::det2x2_capacity(a22, a23, a32, a33)
        );
        m11.assign_det2x2(a22, a23, a32, a33);

        let mut m12 = Expansion::with_capacity(
            Self::det2x2_capacity(a21, a23, a31, a33)
        );
        m12.assign_det2x2(a21, a23, a31, a33);

        let mut m13 = Expansion::with_capacity(
            Self::det2x2_capacity(a21, a22, a31, a32)
        );
        m13.assign_det2x2(a21, a22, a31, a32);

        // 2) form the three products
        let mut t1 = Expansion::with_capacity(
            Self::product_capacity(a11, &m11)
        );
        t1.assign_product(a11, &m11);

        let mut t2 = Expansion::with_capacity(
            Self::product_capacity(a12, &m12)
        );
        t2.assign_product(a12, &m12);

        let mut t3 = Expansion::with_capacity(
            Self::product_capacity(a13, &m13)
        );
        t3.assign_product(a13, &m13);

        // 3) combine: (t1 - t2) + t3
        let mut tmp = Expansion::with_capacity(self.capacity());
        tmp.assign_sum(&t1, &t3);
        self.assign_diff(&tmp, &t2)
    }

    /// Compute the capacity needed to form the 2×2 determinant
    /// of the four expansions a11, a12, a21, a22.
    ///
    /// Mirrors C++:
    /// ```cpp
    /// return product_capacity(a11, a22)
    ///      + product_capacity(a21, a12);
    /// ```
    pub fn det2x2_capacity(
        a11: &Expansion, a12: &Expansion,
        a21: &Expansion, a22: &Expansion
    ) -> usize {
        Self::product_capacity(a11, a22)
            + Self::product_capacity(a21, a12)
    }

    /// Assign to `self` the 2×2 determinant of the four expansions:
    ///     a11·a22 − a12·a21
    ///
    /// # Panics
    /// Panics if `self.capacity()` is less than `det2x2_capacity(...)`.
    pub fn assign_det2x2(
        &mut self,
        a11: &Expansion, a12: &Expansion,
        a21: &Expansion, a22: &Expansion
    ) -> &mut Self {
        // build product a11 * a22
        let mut p1 = Expansion::with_capacity(
            Self::product_capacity(a11, a22)
        );
        p1.assign_product(a11, a22);

        // build product a12 * a21
        let mut p2 = Expansion::with_capacity(
            Self::product_capacity(a12, a21)
        );
        p2.assign_product(a12, a21);

        // self = p1 - p2
        self.assign_diff(&p1, &p2)
    }
}

impl Expansion {
    /// Return the number of components in the expansion.
    pub fn length(&self) -> usize {
        self.data.len()
    }

    /// Compute capacity needed to multiply two expansions a and b.
    pub fn product_capacity(a: &Expansion, b: &Expansion) -> usize {
        a.length().saturating_mul(b.length()).saturating_mul(2)
    }

    /// Assign `self` = a + b (expansion sum).
    /// Naively concatenates the two expansions and then calls `optimize`.
    pub fn assign_sum(&mut self, a: &Expansion, b: &Expansion) -> &mut Self {
        let new_len = a.length() + b.length();
        if self.data.capacity() < new_len {
            self.data.reserve(new_len - self.data.capacity());
        }

        self.data.extend_from_slice(a.data());
        self.data.extend_from_slice(b.data());

        self.data.shrink_to_fit();
        self.optimize();
        self
    }

    /// Assign `self` = a - b (expansion difference).
    pub fn assign_diff(&mut self, a: &Expansion, b: &Expansion) -> &mut Self {
        // build negated b in a temp
        let mut nb = Expansion { data: b.data.clone() };
        nb.negate();
        self.assign_sum(a, &nb)
    }

    /// Assign `self` = a * b (expansion product).
    /// Uses Shewchuk's two_product on each coefficient pair.
    pub fn assign_product(&mut self, a: &Expansion, b: &Expansion) -> &mut Self {
        let cap = Self::product_capacity(a, b);
        assert!(cap <= self.capacity(), "assign_product: capacity too small");
        let mut idx = 0;
        // for each pair (i,j), compute two_product and write low,high
        for &ai in a.data.iter().take(a.length()) {
            for &bi in b.data.iter().take(b.length()) {
                let (low, high) = two_product(ai, bi);
                self.data[idx] = low;
                self.data[idx + 1] = high;
                idx += 2;
            }
        }
        self.data.shrink_to_fit();
        self.optimize();
        self
    }

    /// Remove trailing zero components to maintain a canonical form.
    ///
    /// After this call, `length()` is the smallest index such that
    /// the last component is non-zero, or zero if all components are zero.
    pub fn optimize(&mut self) -> &mut Self {
        while let Some(&last) = self.data.last() {
            if last == 0.0 {
                self.data.pop();
            } else {
                break;
            }
        }
        self
    }
}

/// Two-product: return (low, high) parts of ai*bi
#[inline]
fn two_product(a: f64, b: f64) -> (f64, f64) {
    let x = a * b;

    (x, two_product_tail(a, b, x))
}

#[inline]
fn two_product_tail(a: f64, b: f64, x: f64) -> f64 {
    let (ahi, alo) = split(a);
    let (bhi, blo) = split(b);

    let err1 = x - (ahi * bhi);
    let err2 = err1 - (alo * bhi);
    let err3 = err2 - (ahi * blo);

    (alo * blo) - err3
}

const SPLITTER: f64 = 134_217_729f64;

#[inline]
fn split(a: f64) -> (f64, f64) {
    let c = SPLITTER * a;
    let abig = c - a;
    let ahi = c - abig;
    let alo = a - ahi;

    (ahi, alo)
}