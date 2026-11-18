"""Exact sorted radius iterates for symbolic real inputs."""
import z3


def network(values,choose):
    values=list(values)
    for i,j in ((0,1),(2,3),(0,2),(1,3),(1,2)):
        a,b=values[i],values[j]
        values[i],values[j]=choose(a>=b,a,b),choose(a>=b,b,a)
    return tuple(values)


def symbolic(n):
    q=tuple(z3.Reals('a b c d'));states=[q]
    domain=[q[0]<1,q[0]>=q[1],q[1]>=q[2],q[2]>=q[3],q[3]>=0]
    equations=[]
    for j in range(n):
        x=states[-1]
        values=(x[0],*(z3.If(x[i]+x[i+1]>=1,x[i]+x[i+1]-1,1-x[i]-x[i+1]) for i in range(3)))
        ordered=network(values,z3.If)
        target=tuple(z3.Reals(' '.join(f'r{j+1}_{i}' for i in range(4))))
        equations.extend(a==b for a,b in zip(target,ordered));states.append(target)
    return q,states,domain,equations
