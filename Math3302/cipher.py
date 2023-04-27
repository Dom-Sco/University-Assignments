import collections
import numpy as np
import string

cipher = "DGAEKYMMWWRUDPPEIHISVQIKGZZMMGRVYQPTPMMJKQVRERUZTYZTDTDPQIJREBOSEXCZGFLTYQRBMJJSFBCRRDMNRFKPDRXCWIKVIWMQFIMZBGCQXMKKPVRGZXNSZTIGPVUUTWVCNKMGRVYGRDMOSSVKVPCVLKMGZODZVGQMVZGINQPWIOLVIAAEDZTTRRDMJRKYMVYKYWMMKPWALZJBJQSJBJRVRVYQPJQNXVIUMWAFMBEIXMMCNYWHEIIQZHKYMWPRTSNQZKPVWZEMQIIJIRQPWIOLVIWMQPDWOLVIIIHEVDZVJREVRPCQFIEVANSWVQOLVIWAXYVUASIKPZMIUITWNVZZPFEOWIWFZZXYVLVCJFNKLFKWBVRGPNQPWQMWKWIIGZVAMIXRZYMEXECEKKPZCNVZZPZBMRIIVCIVVRAJRRSTTHVIQQIUWZJQKYMDVKFUWWKFVZWKYMNLRGMJJKYMGIKKMMWFEUTJRKPZVJXIQIDVIISUUQYIRKPVXYVEVWRJYPEIVAOSLKLVVBDIIAZKPXYICGWPRTSCEZINMSDKPZGYRZVGKVZVRUKCMRFWBCIZEAXVZGBDSERTNSXVWMKZRVVAZWMJJKYMVFFMMDHIVEVGYZTYMJYKJRTCCNMFEBCEKDGHSKYMMARJNMITBTZHRELNMTBTTXFWQQICZBOPVJBJRVCWUIEXMNIRTPVFFLBVJFFBVRURPVPWCWIKNYQXLNVZZEIIIIKVUQIEEVIOVFNJZWZUMOLVZZBVRMMVRUNMMIJRKMIUKWOLVDMHSIPWAJZMMGMKKTZFIFBCIIJWAQZEMRLFXIQILGBMCZEOOSXVBVPZMQIKVOKZIUZVBPPVIMPPZVOLRKCIMMVZNECJBMYXXTZMRDQIHVSBZHWFZVFVCQZJZIMGMXZWPWCPMIXVIBVMEVLOLRKBCIPYIYECCJZIESWMRFEBCIZIJVGBJEDXYKPZMIYIIHJZVOLVZZOVFLAZVJGWXOVKAVRUYIYRVMMMXRBMIXYVUJYKZVOLZJAOEKVWAIOZAOIETMJYIJEVWKYMHEIJPXSLEBMCUFEIFPKPZVZMMMAZKPDRRJBCIIZDZVNFCIHKNMIXPDQGIJFNOLVJMVQPWQMWKDWNXMZDDHRELWVFRLDQGIMNWZFVJJKYMDHVEBDXPFNOLZEONWVVUNXFDMOSYRDZFVVVBEZEMYSERUZQFIIWPVIIREWKMMRFFVOSNRZYWVMMIMEXIOWLTPVXZDMDJFLVYSLKNJVTVZOEZEBCEKKPDWSCMVOGCIXIFMMMKIFEIAZKPIIKKTZWNRAOLVTPPVTYGVVURVYXYRBKLZCQKTZIZDTCRBZSWKPDWGRZDWYRVYECJWBIFIODEEREDJVFNOLVRJJZVNMMIUVIYEEUJPVZVLVRUKPVXRCMSEEUMMFRIBCSCFUZARSZVLRDBJFZRAVRUIWBIIZVAEEKKCMCUZZRFWBCIRWWMIJRQYAVIMVPJFLZEURVYFLIQZHRELOLRKBCIURZFJCRBRMCUMMRVJAWIPFVYXYVKCYITPTEIUQIXVIAZGKVLRMKYLTOVJIIHDFCIHJRVYKRKMNAZKPNGRKBZVVUKVXKCMAIVUQIKFEQOARJBCIDRZNLVJIIHKYIOXYVTJACVIYIECQIISVGJRUNINXYVZDZVIIIHKYIOXYVLDWKRVOWRMIBICRQMJIFURLZTPOLVNQIHNRAMYJYQIKNRAOLVJMVEEUBCEKKPZWDRTGFLELGIFWACMMVZNKIFEDRXRNMEZUWAMKRTGEEUJZKZEVDRXKWXVPNINTZGPJPUPWPVEFQNITIQZHRKMMVZSTZZFZKZEJRUVRJKIMXVUCKJIFUVQFEOOLVXZVZVJIOXYVADHVFNOLVTPPVTYXJVTYSZIGJBDPCPWPPZKBGIUVDDPFIQGPTLBTSLIBCVFRBVJVRZAYCDIIECCQIGFRZNIXIMTAZKPVKIVIOMIFVJRYZAGIXRUVRNZBCRFYIOEEUEDXY"

def phi(c):
    x = list(collections.Counter(c).items())
    n = len(c)
    phi = 0
    for i in range(len(x)):
        phi += x[i][1]*(x[i][1]-1)
    phi = phi/(n*(n-1))
    return phi

def friedmanMethod2(c):
    n = len(c)
    phi_0 = 1/26
    phi_L = 0.0667
    phi_T = phi(c)
    m = (n*(phi_L-phi_0))/((n-1)*phi_T-n*phi_0+phi_L)
    return m

print(friedmanMethod2(cipher))

def array2string(arr):
    strng = ""
    for i in range(len(arr)):
        strng += arr[i]
    return strng

def string2array(strng):
    array = np.array([])
    for i in range(len(strng)):
        array = np.append(array,strng[i])
    return array

def cipherSplit(m,c):
    split = []
    cipher = []
    for i in range(len(c)):
        cipher.append(c[i])
    cipher = np.array(cipher)
    add = m - (len(cipher)%m)
    for j in range(add):
        cipher = np.append(cipher,np.array([""]))
    cipher = cipher.reshape(int((len(cipher))/m),m)
    for k in range(m):
        split.append(array2string(cipher[:,k]))
    return split

def friedmanMethod1(c,maxm):
    phimat = np.zeros([maxm,maxm])
    loss = np.tril(0.0667*np.ones([maxm,maxm]))
    for i in range(1,maxm+1):
        splits = cipherSplit(i, c)
        for j in range(len(splits)):
            phimat[i-1,j] += phi(splits[j])
    error = np.abs(phimat-loss)
    error = np.sum(error,axis=1)
    m = np.where(error == np.min(error))[0][0]+1
    return m
            
print(friedmanMethod1(cipher,10))


def letterFrequencies(x):
    freq = {}
 
    for i in x:
        if i in freq:
            freq[i] += 1
        else:
            freq[i] = 1
    total = sum(freq.values(), 0.0)
    freq = {k: v / total for k, v in freq.items()}
    return freq
    
def shiftCalculation(c,m):
    d = {chr(i+65):i for i in range(0,26)}
    letters = ["E","T","A","O","I","N","S","H","R","D","L","U","C","M","W","F","G","Y","P","B","V","K","J","Z","X","Q"] #Letters from highest frequency in english language to lowest frequency
    split = cipherSplit(m, c)
    shifts = []
    for i in range(len(split)):
        listletters = []
        freq = letterFrequencies(split[i])
        sort = sorted(freq.values(),reverse=True)
        for j in range(len(freq)):
            listletters.append(d[list(freq.keys())[list(freq.values()).index(sort[j])]])
        distance = []
        for k in range(2):
            ma = max(d[letters[0]],listletters[k])
            mi = min(d[letters[0]],listletters[k])
            distance.append((ma-mi))
        shifts.append(distance)
    return shifts

shifts = shiftCalculation(cipher,5)

#the key is [17,8,21,4,17]

def decrypt(c,m):
    d = {chr(i+65):i for i in range(0,26)}
    split = cipherSplit(m, c)
    shifts = [17,8,21,4,17]
    arraylen = max(len(i) for i in split)
    decrypt = []
    for i in range(len(split)):
        strng = ""
        for j in range(len(split[i])):
            s = (d[split[i][j]]-shifts[i])%26
            strng += list(d.keys())[list(d.values()).index(s)]
        decrypt.append(strng)
    message = string2array(decrypt[0])
    if len(message)!=arraylen:
        while len(message)!=arraylen:
            message = np.append(message,0)
    for k in range(len(decrypt)-1):
        arr = string2array(decrypt[k+1])
        if len(arr)!=arraylen:
            while len(arr)!=arraylen:
                arr = np.append(arr,0)
        message = np.vstack((message,string2array(decrypt[k+1])))
    message = message.T
    shape = message.shape[0]*message.shape[1]
    message = message.reshape(shape)
    return message


message = decrypt(cipher,5)
print(array2string(message))