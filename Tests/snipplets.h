/*
 * snipplets.h
 *
 *  Created on: 18.02.2011
 *      Author: Martin Rueckl
 */

#ifndef SNIPPLETS_H_
#define SNIPPLETS_H_

class A{
	public:
		void bla(){
			cout<<"class A::bla"<<endl;
		}
};

class B{
	public:
		void bla(){
			cout<<"class B::bla"<<endl;
		}
};


template<typename T>
class tester{
	public:
		tester(T _t):t(_t){
			t.bla();
		}
	private:
		T t;
};




void templates(){
	A a;
	B b;
	tester<A> ta(a);
	tester<B> tb(b);
}

class A{
	public:
	void blub(){cout<<"A::blub"<<endl;}
};

class B{
	public:
	void bla(){cout<<"B::bla"<<endl;}
};

class C:public A{
	public:
		void blub(){cout <<"C::blub"<<endl; }
};

void static_cast_kram()
{
	A a;
	//B b = static_cast<B>(a);
	//b.bla();
	C c;
	A a2 = static_cast<A>(c);
	c.blub();
}

#endif /* SNIPPLETS_H_ */
