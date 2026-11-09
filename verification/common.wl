(* Shared helpers for the Wolfram verification scripts of TOVBrachistochrone v2.

   Every script records its checks with tovRecord and ends with tovFinish. The runner
   (tools/run_verification.py) passes the output path in the environment variable
   TOV_RESULT_JSON and a unique run identifier in TOV_RUN_ID, which is copied into the result so that
   a result written by another run is rejected. The result is written to a temporary file and renamed,
   so an interrupted run never leaves a half-written JSON.

   Evidence kinds:
     "exact"             symbolic identity or decision (Simplify, Reduce, SolveAlways, polynomial algebra)
     "numerical_check"   a quantitative check of a stated result at given parameters, with a stated tolerance:
                         a closed form against an independent quadrature or ODE, or the evaluation of an
                         analytic statement (a sign, a bracket, an asymptotic ratio) at sample points
     "diagnostic"        a sampled numerical statement that the paper labels as diagnostic, not as proof
     "negative_control"  a deliberate mutation; status "pass" means the mutation was detected

   Exit status: 0 if every check passes, 1 otherwise, 2 if the script aborts on a violated assumption. *)

tovSchemaVersion = "tov-result/1";
tovChecks = {};
tovValues = <||>;

tovRecord[id_String, kind_String, description_String, ok_] := Module[{status},
   If[! MemberQ[{"exact", "numerical_check", "diagnostic", "negative_control"}, kind],
    Print["ABORT: unknown evidence kind ", kind]; Exit[2]];
   If[MemberQ[#["id"] & /@ tovChecks, id], Print["ABORT: duplicate check id ", id]; Exit[2]];
   status = If[TrueQ[ok], "pass", "fail"];
   AppendTo[tovChecks, <|"id" -> id, "kind" -> kind, "description" -> description, "status" -> status|>];
   Print[If[status == "pass", "PASS  ", "FAIL  "], id, "  (", kind, ")  ", description]];

(* numbers are stored as decimal strings ("0.123", "1.5e-13") so that high-precision values survive JSON;
   the Wolfram exponent mark "*^" is rewritten as "e" *)
tovValue[key_String, value_?NumericQ, digits_Integer: 30] :=
  (tovValues[key] = StringReplace[ToString[N[value, digits], InputForm, NumberMarks -> False], "*^" -> "e"]);
tovValue[key_String, value_String] := (tovValues[key] = value);
tovValue[key_String, value_] := (tovValues[key] = ToString[value, InputForm]);

(* canonical identifier of an exact rational parameter: 1/10 -> "1_10", 8888/10000 -> "1111_1250", 2 -> "2" *)
tovTag[r_Rational] := ToString[Numerator[r]] <> "_" <> ToString[Denominator[r]];
tovTag[r_Integer] := ToString[r];
tovTag[r_] := (Print["ABORT: tovTag needs an exact rational parameter, got ", r]; Exit[2]);

tovAbort[message_String] := (Print["ABORT: ", message]; Exit[2]);

tovFinish[jobId_String] := Module[{out = Environment["TOV_RESULT_JSON"], passed, total, result, tmp},
   total = Length[tovChecks];
   passed = Count[tovChecks, _?(#["status"] == "pass" &)];
   result = <|"schema_version" -> tovSchemaVersion, "job_id" -> jobId,
     "run_id" -> Replace[Environment["TOV_RUN_ID"], $Failed -> "unset"],
     "engine" -> <|"name" -> "Wolfram Language", "version" -> $Version|>,
     "checks" -> tovChecks, "values" -> tovValues,
     "counts" -> <|"total" -> total, "passed" -> passed, "failed" -> total - passed|>|>;
   If[StringQ[out],
    tmp = out <> ".tmp";
    Export[tmp, result, "JSON"];
    RenameFile[tmp, out, OverwriteTarget -> True]];
   Print["passed ", passed, " of ", total];
   Exit[If[passed == total && total > 0, 0, 1]]];
