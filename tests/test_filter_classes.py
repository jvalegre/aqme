from datetime import date
import unittest
import random
from datetime import datetime
from pyconfort.filter import Filter, CompoundFilter

random.seed(datetime.now())

def f1(x):
    return x > 5.0
def f2(x): 
    return len(x) <= 6
def f3(x): 
    return x < 15.0 


class TestFilter(unittest.TestCase):
    @classmethod
    def setUpClass(cls): 
        cls.dataset1 = list(range(10))
        cls.dataset2 = ['dog','cat','seal','eagle','rabbit',
                        'narwhal','elephant','groundhog']
        cls.dataset3 = [(7,'apple'),
                        (2,'pears'),
                        (9,'bananas'),
                        (4,'pineapples'),
                        (6,'fruits')]
    def test_init(self):
        msg = 'Wrong initialization'
        A = Filter()
        self.assertTrue(A.function('dummy'),msg)
        self.assertIsNone(A.dataset,msg)
        self.assertIsNone(A.outcomes,msg)
        self.assertIsNone(A._discarded,msg)
        self.assertIsNone(A._accepted,msg)
        f = lambda x: False
        B = Filter(f)
        self.assertIsNotNone(B.function,msg)
        self.assertFalse(B.function('dummy'),msg)
    def test_calc_outcomes(self):
        A = Filter(f1) 
        B = Filter(f2)
        test1 = A.calc_outcomes(self.dataset1,key=lambda x: x)
        test2 = B.calc_outcomes(self.dataset2,key=lambda x: x)
        test3 = A.calc_outcomes(self.dataset3,key=lambda x: x[0])
        out1 = (False,False,False,False,False,False,True,True,True,True)
        out2 = (True,True,True,True,True,False,False,False)
        out3 = (True,False,True,False,True)
        self.assertEqual(test1,out1)
        self.assertEqual(test2,out2)
        self.assertEqual(test3,out3)
    def test_apply(self):
        A = Filter(f1) 
        B = Filter(f2)
        dataset1 = list(self.dataset1)
        dataset2 = list(self.dataset2)
        dataset3 = list(self.dataset3)

        random.shuffle(dataset1)
        random.shuffle(dataset2)
        random.shuffle(dataset3)

        A.apply(dataset1)
        self.assertEqual(dataset1,A.dataset)
        self.assertIsNotNone(A.outcomes)
        B.apply(dataset2)
        self.assertEqual(dataset2,B.dataset)
        self.assertIsNotNone(B.outcomes)

        old_outcomes = A.outcomes
        key = lambda x: x[0]
        with self.assertRaises(ValueError): 
            A.apply(dataset3,key=key)
        
        A.apply(dataset3,key=key,force=True)
        self.assertNotEqual(set(dataset1),set(A.dataset))
        self.assertEqual(set(dataset3),set(A.dataset))
        self.assertNotEqual(A.outcomes,old_outcomes)

    def test_extend(self):
        dataset1 = list(self.dataset1)
        dataset3 = list(self.dataset3)
        accepted = [6,7,8,9,(7,'apple'),(9,'bananas'),(6,'fruits')]
        random.shuffle(dataset1) 
        A = Filter(f1)
        A.apply(dataset1)
        A.extend(dataset3,key=lambda x: x[0])
        self.assertEqual(A.dataset[:10],dataset1)
        self.assertEqual(A.dataset[10:],dataset3)
        self.assertEqual(set(A.accepted[:4]),set(accepted[:4]))
        self.assertEqual(A.accepted[4:],accepted[4:])

    def test_clean(self):
        dataset1 = list(self.dataset1)
        A = Filter(f1)
        self.assertIsNone(A.dataset)
        self.assertIsNone(A.outcomes)
        A.apply(dataset1)
        self.assertIsNotNone(A.dataset)
        self.assertIsNotNone(A.outcomes)
        A.clean()
        self.assertIsNone(A.dataset)
        self.assertIsNone(A.outcomes)
    def test_properties(self):
        accepted = [(7,'apple'),
                    (2,'pears'),
                    (6,'fruits')]
        discarded = [(9,'bananas'),
                     (4,'pineapples')]
        dataset3 = list(self.dataset3)
        random.shuffle(dataset3)
        A = Filter(f2)
        self.assertIsNone(A.accepted)
        self.assertIsNone(A.discarded)
        A.apply(dataset3,key=lambda x: x[1])
        A_accepted = A.accepted
        A_discarded = A.discarded
        self.assertEqual(set(A_accepted),set(accepted))
        self.assertEqual(set(A_discarded),set(discarded))
        # Check the cache works
        A.dataset.append('Fake')
        A.outcomes = A.outcomes + (True,)
        self.assertEqual(A.accepted,A_accepted)
        self.assertEqual(A.discarded,A_discarded)
        # Check the cleaning resets them 
        A.clean()
        self.assertIsNone(A.accepted)
        self.assertIsNone(A.discarded)

class TestCompoundFilter(unittest.TestCase):
    @classmethod
    def setUpClass(cls): 
        cls.dataset1 = list(range(20))
    def test_init(self):
        msg = 'Wrong initialization'
        A = Filter()
        f = lambda x: False
        B = Filter(f)
        C = CompoundFilter(A)
        D = CompoundFilter(A,B) 
        E = CompoundFilter()
        self.assertTrue(len(C.filters)==1,msg)
        self.assertTrue(len(D.filters)==2,msg)
        self.assertTrue(len(E.filters)==0,msg)
        self.assertIsNone(D.dataset,msg)
        self.assertIsNone(D.outcomes,msg)
        self.assertIsNone(D._discarded,msg)
        self.assertIsNone(D._accepted,msg)
        self.assertTrue(D.function('dummy'),msg)
        self.assertTrue(E.function(None),msg)
    def test_calc_outcomes(self):
        A = Filter(f1) 
        B = Filter(f3)
        C = CompoundFilter(A,B)
        test = C.calc_outcomes(self.dataset1,key=lambda x: x)
        outA = tuple(i=='1' for i in '00000011111111111111')
        outB = tuple(i=='1' for i in '11111111100000')
        self.assertEqual(test,[outA,outB])
    def test_apply(self):
        A = Filter(f1) 
        B = Filter(f3)
        C = CompoundFilter(A,B)
        dataset1 = list(self.dataset1)
        dataset2 = list(self.dataset1)
        
        random.shuffle(dataset1)
        random.shuffle(dataset2)

        C.apply(dataset1)
        self.assertSetEqual(set(dataset1),set(C.dataset))
        self.assertNotEqual(set(dataset1),set(B.dataset))
        self.assertIsNotNone(C.outcomes)

        old_outcomes = C.outcomes
        with self.assertRaises(ValueError): 
            C.apply(dataset2)
        
        C.apply(dataset2,force=True)
        self.assertNotEqual(dataset1,C.dataset)
        self.assertEqual(dataset2,C.dataset)
        self.assertNotEqual(C.outcomes,old_outcomes)
    def test_extend(self):
        A = Filter(f1) 
        B = Filter(f3)
        C = CompoundFilter(A,B)
        dataset1 = list(self.dataset1[:10])
        dataset2 = list(self.dataset1[10:])
        accepted = [i for i,t in enumerate('000000111111111000000') if t=='1']        
        
        random.shuffle(dataset1)

        C.apply(dataset1)
        C.extend(dataset2)
        self.assertEqual(C.dataset,dataset1+dataset2)
        self.assertEqual(C.dataset[:10],dataset1)
        self.assertEqual(set(C.accepted[:4]),set(accepted[:4]))
        self.assertEqual(C.accepted[4:],accepted[4:])
    def test_clean(self):
        dataset1 = list(self.dataset1)
        A = Filter(f1)
        B = Filter(f3)
        C = CompoundFilter(A,B)
        self.assertIsNone(A.dataset)
        self.assertIsNone(A.outcomes)
        self.assertIsNone(B.dataset)
        self.assertIsNone(B.outcomes)
        self.assertIsNone(C.dataset)
        self.assertIsNone(C.outcomes)
        C.apply(dataset1)
        self.assertIsNotNone(A.dataset)
        self.assertIsNotNone(A.outcomes)
        self.assertIsNotNone(B.dataset)
        self.assertIsNotNone(B.outcomes)
        self.assertIsNotNone(C.dataset)
        self.assertIsNotNone(C.outcomes)
        C.clean()
        self.assertIsNone(A.dataset)
        self.assertIsNone(A.outcomes)
        self.assertIsNone(B.dataset)
        self.assertIsNone(B.outcomes)
        self.assertIsNone(C.dataset)
        self.assertIsNone(C.outcomes)
    def test_properties(self):
        A = Filter(f1) 
        B = Filter(f3)
        C = CompoundFilter(A,B)
        outcome = '00000011111111100000'
        accepted = tuple(i for i,t in enumerate(outcome) if t=='1')
        discarded = tuple(i for i,t in enumerate(outcome) if t=='0')
        dataset1 = list(self.dataset1)
        
        random.shuffle(dataset1)

        self.assertIsNone(A.accepted)
        self.assertIsNone(A.discarded)
        self.assertIsNone(B.accepted)
        self.assertIsNone(B.discarded)
        self.assertIsNone(C.accepted)
        self.assertIsNone(C.discarded)
        C.apply(dataset1)
        A_accepted = set(A.accepted)
        B_accepted = set(B.accepted)
        C_accepted = set(C.accepted)
        A_discarded = set(A.discarded)
        B_discarded = set(B.discarded)
        C_discarded = set(C.discarded)
        self.assertEqual(C_accepted,A_accepted and (B_accepted))
        self.assertEqual(C_accepted,set(accepted))
        self.assertEqual(C_discarded,A_discarded | B_discarded)
        self.assertEqual(C_discarded,set(discarded))
        # Check the cache works
        C.outcomes.append((True,)*len(C.outcomes[-1]))
        self.assertEqual(set(C.accepted),C_accepted)
        self.assertEqual(set(C.discarded),C_discarded)
        # Check the cleaning resets them 
        C.clean()
        self.assertIsNone(A.accepted)
        self.assertIsNone(A.discarded)
        self.assertIsNone(B.accepted)
        self.assertIsNone(B.discarded)
        self.assertIsNone(C.accepted)
        self.assertIsNone(C.discarded)

    def test_accepted_from(self):
        A = Filter(f1) 
        B = Filter(f3)
        C = CompoundFilter(A,B)
        dataset1 = list(self.dataset1)
        accepted = [[i for i,t in enumerate('00000011111111111111') if t=='1'],
                    [i for i,t in enumerate('00000011111111100000') if t=='1']]
        C.apply(dataset1)
        
        self.assertEqual(C.accepted_from(0),accepted[0])
        self.assertEqual(C.accepted_from(1),accepted[1])
    def test_discarded_from(self):
        A = Filter(f1) 
        B = Filter(f3)
        C = CompoundFilter(A,B)
        dataset1 = list(self.dataset1)
        discarded= [[i for i,t in enumerate('00000011111111111111') if t=='0'],
                    [i for i,t in enumerate('00000011111111100000') if t=='0']]
        C.apply(dataset1)
        
        self.assertEqual(C.discarded_from(0),discarded[0])
        self.assertEqual(C.discarded_from(1),discarded[1])
