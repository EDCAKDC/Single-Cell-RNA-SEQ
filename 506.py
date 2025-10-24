class Solution:
    def checkPerfectNumber(self, num: int) -> bool:
        if num <= 1:
            return False
        n=1
        m=2
        while m*m < num:
            if num%m==0:
                n += m
                l = num //m
                if l != m:
                    n += l
            m += 1
        return n==num