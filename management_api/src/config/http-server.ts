import { api } from "../sdkgen/api-generated";
import { SdkgenHttpServer } from "@sdkgen/node-runtime";
import "../http/server";
import { QUERY } from "../sql";

export function runHttpServer() {
	// run the server
	const server = new SdkgenHttpServer(api, {});
	server.listen(8000);

	api.hook.onRequestStart = async function (ctx) {
		if (ctx.request.name !== "auth" && ctx.request.name !== "createAccount") {
			const userId = await QUERY.users.getUserIdByDeviceId(
				ctx.request.deviceInfo.id,
			);
			if (!userId) throw new Error("Não autenticado");
			ctx.request.extra["userId"] = userId;
			return null;
		}
		return null;
	};
}
