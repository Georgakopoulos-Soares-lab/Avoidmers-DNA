import re
aba_pattern = re.compile(r'(.+)(.+)\1')
sequences = set()


def is_aba(x):
    N = len(x)
    if N <= 2:
        return False
        
    return any(x[:i] == x[-i:] for i in range(1, (len(x)+1)//2))
    
def zimin_density(x):
    total_aba = 0 
    N = len(x)
    total_s = N * (N+1) // 2
    
    for l in range(1, N+1):
        for i in range(N-l+1):
            chunk = x[i:i+l]
            
            # is_aba = re.search(aba_pattern, chunk)
            if is_aba(chunk):
                total_aba += 1
                
    return total_aba, total_s, total_aba / total_s


if __name__ == "__main__":


    pass
