GiD Post Results File 1.2
# meshio++: a hand-authored fixture, not writer output
GaussPoints "gp_surf" ElemType Triangle "surf"
Number Of Gauss Points: 1
Natural Coordinates: Internal
End GaussPoints
Result "T" "analysis" 1 Scalar OnNodes
Values
1 10
2 20
3 30
4 40
End Values
Result "T" "analysis" 2 Scalar OnNodes
Values
1 100
2 200
3 300
4 400
End Values
Result "flux" "analysis" 1 Vector OnNodes
Values
1 1 2 3
2 4 5 6
3 7 8 9
4 10 11 12
End Values
Result "q" "analysis" 1 Scalar OnGaussPoints "gp_surf"
Values
1 5
2 6
End Values
GaussPoints "gp2_surf" ElemType Triangle "surf"
Number Of Gauss Points: 2
Natural Coordinates: Internal
End GaussPoints
Result "twogp" "analysis" 1 Scalar OnGaussPoints "gp2_surf"
Values
1 5
 6
2 7
 8
End Values
