Print("GAP_VERSION=", GAPInfo.Version, "\n");

F := FreeGroup("g1", "g2");;
g1 := F.1;;
g2 := F.2;;
G := F / [g1^-4 * g2 * g1 * g2,
          g1 * g2^-2 * g1 * g2];;

Print("SIZE=", Size(G), "\n");
Print("PERFECT=", IsPerfectGroup(G), "\n");
Print("ABELIAN_INVARIANTS=", AbelianInvariants(G), "\n");
QUIT;
