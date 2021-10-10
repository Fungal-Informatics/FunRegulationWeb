import { LOGIC } from "../logic";
import { api } from "../sdkgen/api-generated";

api.fn.auth = async (ctx, { email, password }) => {
	await LOGIC.users.authUser(ctx.request.deviceInfo.id, email, password);
};

api.fn.logout = async (ctx) => {
	await LOGIC.users.deauthUser(ctx.request.deviceInfo.id);
};

api.fn.myUser = async (ctx) => {
	return LOGIC.users.getUser(ctx.request.extra.userId);
};

api.fn.createAccount = async (_ctx, { newAcount }) => {
	await LOGIC.users.createUser(
		newAcount.name,
		newAcount.email,
		newAcount.password,
	);
};
