/*
 * Downloaded from Heng Li's github repository 
 * URL : https://gist.github.com/lh3/974ced188be2f90422cc
 */

/*
  For any 1<k<=64, let mask=(1<<k)-1. hash_64() is a bijection on [0,1<<k), which means
  hash_64(x, mask)==hash_64(y, mask) if and only if x==y. hash_64i() is the inversion of
  hash_64(): hash_64i(hash_64(x, mask), mask) == hash_64(hash_64i(x, mask), mask) == x.
*/

namespace conn
{
  namespace graphGen
  {

    // Thomas Wang's integer hash functions. See <https://gist.github.com/lh3/59882d6b96166dfc3d8d> for a snapshot.
    // NOTE: T should be a 64 bit type
    template <typename T>
      void hash_64(T &key)
      {
        key = (~key + (key << 21)); // key = (key << 21) - key - 1;
        key = key ^ key >> 24;
        key = ((key + (key << 3)) + (key << 8)); // key * 265
        key = key ^ key >> 14;
        key = ((key + (key << 2)) + (key << 4)); // key * 21
        key = key ^ key >> 28;
        key = (key + (key << 31));
      }

    // The inversion of hash_64(). Modified from <https://naml.us/blog/tag/invertible>
    template <typename T>
      void hash_64i(T &key)
      {
        uint64_t tmp;

        // Invert key = key + (key << 31)
        tmp = (key - (key << 31));
        key = (key - (tmp << 31));

        // Invert key = key ^ (key >> 28)
        tmp = key ^ key >> 28;
        key = key ^ tmp >> 28;

        // Invert key *= 21
        key = (key * 14933078535860113213ull);

        // Invert key = key ^ (key >> 14)
        tmp = key ^ key >> 14;
        tmp = key ^ tmp >> 14;
        tmp = key ^ tmp >> 14;
        key = key ^ tmp >> 14;

        // Invert key *= 265
        key = (key * 15244667743933553977ull);

        // Invert key = key ^ (key >> 24)
        tmp = key ^ key >> 24;
        key = key ^ tmp >> 24;

        // Invert key = (~key) + (key << 21)
        tmp = ~key;
        tmp = ~(key - (tmp << 21));
        tmp = ~(key - (tmp << 21));
        key = ~(key - (tmp << 21));
      }
  }
}
